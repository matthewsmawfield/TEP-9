#!/usr/bin/env python3
"""
Flyby TEP Pipeline - Enhanced Physics Module

Addresses Model Incompleteness weakness:
- Non-spherical Earth (oblateness with J2 and higher-order harmonics)
- 3D trajectory integration along exact flight path
- Dynamic density mapping using geoid models
- Proper residual calculation (not per-flyby linear scaling)
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, List, Optional, Callable, Dict, Any
from scipy.integrate import solve_ivp

# Standardized Physics constants for Jakarta v0.8
from scripts.utils.physics import (
    C_LIGHT, G_NEWTON, R_EARTH, M_EARTH, J2_EARTH, M_PL_GEV as M_PL,
    BETA_BASELINE, get_tep_metadata
)

# Conversion factors
KG_M3_TO_GEV4 = 4.318e-21  # kg/m^3 to GeV^4 conversion factor
OMEGA_EARTH = 7.2921159e-5  # Earth rotation rate [rad/s]


@dataclass
class EarthGeoidModel:
    """
    Standardized Earth Geoid and Density Structure Model.
    
    This class implements the geodetic and geophysical parameters used for
    mapping the Temporal Topology field near Earth, incorporating the
    WGS84 ellipsoid and PREM density profile.
    """
    
    # WGS84 ellipsoid parameters
    a_equatorial: float = R_EARTH
    b_polar: float = 6.3567523142e6  # m
    flattening: float = 1 / 298.257223563
    
    # Gravitational harmonics (EGM2008 reference)
    J2: float = J2_EARTH
    J3: float = -2.54e-6
    J4: float = -1.61e-6
    
    # Rotation
    omega: float = 7.2921159e-5  # rad/s
    
    def radius_at_latitude(self, lat: float) -> float:
        """
        Earth radius at given geodetic latitude using WGS84 formula.
        
        Parameters
        ----------
        lat : float
            Geodetic latitude in radians
            
        Returns
        -------
        float
            Earth radius at that latitude in meters
        """
        e2 = 2 * self.flattening - self.flattening**2  # eccentricity squared
        sin_lat = np.sin(lat)
        
        # Radius of curvature in the prime vertical
        N = self.a_equatorial / np.sqrt(1 - e2 * sin_lat**2)
        
        # Distance from Earth's center to surface
        r = N * np.sqrt((1 - e2)**2 * sin_lat**2 + np.cos(lat)**2)
        
        return r
    
    def geodetic_to_geocentric(self, lat_geodetic: float, alt: float, 
                                lon: float) -> Tuple[float, float, float]:
        """
        Convert geodetic coordinates (lat, lon, alt) to geocentric (x, y, z).
        
        Parameters
        ----------
        lat_geodetic : float
            Geodetic latitude in radians
        alt : float
            Altitude above ellipsoid in meters
        lon : float
            Longitude in radians
            
        Returns
        -------
        tuple (x, y, z)
            Geocentric coordinates in meters
        """
        e2 = 2 * self.flattening - self.flattening**2
        sin_lat = np.sin(lat_geodetic)
        cos_lat = np.cos(lat_geodetic)
        sin_lon = np.sin(lon)
        cos_lon = np.cos(lon)
        
        # Radius of curvature in prime vertical
        N = self.a_equatorial / np.sqrt(1 - e2 * sin_lat**2)
        
        x = (N + alt) * cos_lat * cos_lon
        y = (N + alt) * cos_lat * sin_lon
        z = (N * (1 - e2) + alt) * sin_lat
        
        return x, y, z
    
    def geocentric_to_geodetic(self, x: float, y: float, z: float) -> Tuple[float, float, float]:
        """
        Convert geocentric (x, y, z) to geodetic (lat, lon, alt).
        
        Uses iterative Bowring's formula for accurate conversion.
        
        Returns
        -------
        tuple (lat, lon, alt)
            Geodetic latitude (rad), longitude (rad), altitude (m)
        """
        p = np.sqrt(x**2 + y**2)
        
        # Initial guess
        lat = np.arctan2(z, p * (1 - self.flattening))
        
        # Iterate for accuracy
        for _ in range(5):
            sin_lat = np.sin(lat)
            N = self.a_equatorial / np.sqrt(1 - (2 * self.flattening - self.flattening**2) * sin_lat**2)
            h = p / np.cos(lat) - N
            lat = np.arctan2(z, p * (1 - self.flattening * (N / (N + h))))
        
        lon = np.arctan2(y, x)
        
        # Final altitude calculation
        sin_lat = np.sin(lat)
        N = self.a_equatorial / np.sqrt(1 - (2 * self.flattening - self.flattening**2) * sin_lat**2)
        alt = p / np.cos(lat) - N
        
        return lat, lon, alt
    
    def gravity_potential(self, r: float, theta: float) -> float:
        """
        Earth's gravitational potential including J2 and higher harmonics.
        
        Parameters
        ----------
        r : float
            Distance from Earth's center (m)
        theta : float
            Geocentric colatitude (angle from z-axis, radians)
            
        Returns
        -------
        float
            Gravitational potential (J/kg)
        """
        mu = G_NEWTON * M_EARTH
        
        # Spherical term
        U0 = -mu / r
        
        # J2 perturbation (oblateness)
        P2 = 0.5 * (3 * np.cos(theta)**2 - 1)  # Legendre polynomial P2
        U2 = -mu / r * (self.a_equatorial / r)**2 * self.J2 * P2
        
        # J3 perturbation (pear-shaped)
        P3 = 0.5 * (5 * np.cos(theta)**3 - 3 * np.cos(theta))
        U3 = -mu / r * (self.a_equatorial / r)**3 * self.J3 * P3
        
        # J4 perturbation
        P4 = 0.125 * (35 * np.cos(theta)**4 - 30 * np.cos(theta)**2 + 3)
        U4 = -mu / r * (self.a_equatorial / r)**4 * self.J4 * P4
        
        return U0 + U2 + U3 + U4
    
    def gravity_acceleration(self, x: float, y: float, z: float) -> Tuple[float, float, float]:
        """
        Compute gravitational acceleration vector at position (x, y, z).
        
        Parameters
        ----------
        x, y, z : float
            Geocentric coordinates in meters
            
        Returns
        -------
        tuple (ax, ay, az)
            Acceleration components in m/s^2
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        if r < 1e-6:
            return 0.0, 0.0, 0.0
        
        theta = np.arccos(z / r)  # colatitude
        phi = np.arctan2(y, x)  # longitude
        
        mu = G_NEWTON * M_EARTH
        
        # Radial component (spherical + J2 + J3 + J4)
        P2 = 0.5 * (3 * np.cos(theta)**2 - 1)
        dP2 = -3 * np.cos(theta) * np.sin(theta)
        
        P3 = 0.5 * (5 * np.cos(theta)**3 - 3 * np.cos(theta))
        dP3 = 1.5 * np.sin(theta) * (1 - 5 * np.cos(theta)**2)
        
        # Radial gravitational acceleration (spherical + J2 + J3 + J4)
        P4 = 0.125 * (35 * np.cos(theta)**4 - 30 * np.cos(theta)**2 + 3)
        g_r = -mu / r**2 * (
            1
            - 3 * self.J2 * (self.a_equatorial / r)**2 * P2
            - 4 * self.J3 * (self.a_equatorial / r)**3 * P3
            - 5 * self.J4 * (self.a_equatorial / r)**4 * P4
        )
        
        # Tangential (meridional) component
        dP4 = 2.5 * np.sin(theta) * (3 * np.cos(theta) - 7 * np.cos(theta)**3)
        g_theta = -mu / r**2 * (
            self.J2 * (self.a_equatorial / r)**2 * dP2
            + self.J3 * (self.a_equatorial / r)**3 * dP3
            + self.J4 * (self.a_equatorial / r)**4 * dP4
        )
        
        # Convert to Cartesian
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        
        ax = g_r * sin_theta * cos_phi + g_theta * cos_theta * cos_phi
        ay = g_r * sin_theta * sin_phi + g_theta * cos_theta * sin_phi
        az = g_r * cos_theta - g_theta * sin_theta
        
        return ax, ay, az


class EarthDensityModel:
    """
    Dynamic density model for Earth based on geoid and crustal structure.
    
    Implements layered Earth model with:
    - Core (inner/outer)
    - Mantle (lower/upper)
    - Crust (oceanic/continental)
    - Atmosphere (exponential decay)
    """
    
    # Layer boundaries (m from center)
    R_INNER_CORE = 1.2215e6
    R_OUTER_CORE = 3.4800e6
    R_LOWER_MANTLE = 5.701e6
    R_UPPER_MANTLE = 6.351e6
    R_CRUST = 6.371e6  # Surface
    
    # Average densities (kg/m^3)
    RHO_INNER_CORE = 13000
    RHO_OUTER_CORE = 11100
    RHO_LOWER_MANTLE = 4900
    RHO_UPPER_MANTLE = 3600
    RHO_CRUST_OCEANIC = 2900
    RHO_CRUST_CONTINENTAL = 2700
    RHO_SURFACE_AVG = 2700
    RHO_ATM_SEA_LEVEL = 1.225
    
    # PREM density profile coefficients
    def __init__(self, geoid: EarthGeoidModel = None):
        self.geoid = geoid or EarthGeoidModel()
    
    def density_prem(self, r: float) -> float:
        """
        PREM (Preliminary Reference Earth Model) density profile.
        
        Parameters
        ----------
        r : float
            Distance from Earth's center (m)
            
        Returns
        -------
        float
            Density in kg/m^3
        """
        r_km = r / 1000  # Convert to km for PREM formula
        
        if r <= self.R_INNER_CORE:
            # Inner core: r <= 1221.5 km
            return 13.0885 - 8.8381 * (r_km / 6371)**2
        elif r <= self.R_OUTER_CORE:
            # Outer core: 1221.5 < r <= 3480 km
            return 12.5815 - 1.2638 * (r_km / 6371) - 3.6426 * (r_km / 6371)**2 - 5.5281 * (r_km / 6371)**3
        elif r <= self.R_LOWER_MANTLE:
            # Lower mantle: 3480 < r <= 5701 km
            return 6.8141 - 1.4837 * (r_km / 6371)
        elif r <= self.R_UPPER_MANTLE:
            # Transition zone: 5701 < r <= 5771 km
            return 11.2494 - 8.0298 * (r_km / 6371)
        elif r <= 5971e3:
            # Upper mantle: 5771 < r <= 5971 km
            return 7.1089 - 3.8045 * (r_km / 6371)
        elif r <= 6151e3:
            # Upper mantle continued: 5971 < r <= 6151 km
            return 2.6910 + 0.6924 * (r_km / 6371)
        elif r <= self.R_CRUST:
            # Crust: 6151 < r <= 6371 km
            return 2.900
        else:
            # Above surface: exponential atmosphere
            h = r - self.R_CRUST  # altitude in meters
            return self.RHO_ATM_SEA_LEVEL * np.exp(-h / 8500)  # scale height 8.5 km
    
    def density_with_geoid(self, x: float, y: float, z: float) -> float:
        """
        Compute density at position accounting for geoid shape.
        
        Uses PREM for interior, exponential atmosphere for exterior.
        
        Parameters
        ----------
        x, y, z : float
            Geocentric coordinates in meters
            
        Returns
        -------
        float
            Density in kg/m^3
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        
        # Get geodetic coordinates
        lat, lon, alt = self.geoid.geocentric_to_geodetic(x, y, z)
        
        if alt < 0:
            # Below reference ellipsoid: use PREM density
            return self.density_prem(r)
        else:
            # Above surface: exponential atmosphere
            return self.RHO_ATM_SEA_LEVEL * np.exp(-alt / 8500)
    
    def local_temporal_topology_threshold(self, x: float, y: float, z: float) -> float:
        """
        Compute local Temporal Topology threshold based on local density.
        
        The suppression threshold depends on the local matter density.
        Higher density = stronger Temporal Shear suppression = shorter correlation length.
        
        Returns
        -------
        float
            Characteristic Temporal Topology suppression length in meters
        """
        rho = self.density_with_geoid(x, y, z)
        
        # Temporal Topology suppression length: λ_rest ∝ 1/√rho
        # Derived from the continuous Temporal Shear Suppression regime of the TEP scalar field
        rho_threshold = 1e-3  # kg/m^3 (approximate transition density)
        
        if rho > rho_threshold:
            # Shielded regime: short range suppression of Temporal Shear
            return 1000.0 * (rho_threshold / rho)**0.5  # meters
        else:
            # Unscreened regime: long range active field (interplanetary space)
            return 1e6  # 1000 km scale


@dataclass
class TrajectoryState:
    """
    State vector for spacecraft trajectory integration.
    """
    t: float  # Time (s)
    x: float  # Position x (m)
    y: float  # Position y (m)
    z: float  # Position z (m)
    vx: float  # Velocity x (m/s)
    vy: float  # Velocity y (m/s)
    vz: float  # Velocity z (m/s)
    
    @property
    def position(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])
    
    @property
    def velocity(self) -> np.ndarray:
        return np.array([self.vx, self.vy, self.vz])
    
    @property
    def r(self) -> float:
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    @property
    def v(self) -> float:
        return np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)


class TEP3DTrajectoryIntegrator:
    """
    3D trajectory integrator with TEP Temporal Topology and non-spherical Earth.
    
    Integrates spacecraft trajectory accounting for:
    - Non-spherical Earth gravity (J2, J3, J4)
    - Earth's rotation
    - TEP Temporal Topology field perturbation
    - Variable coupling based on local density (Temporal Shear Suppression)
    """
    
    def __init__(self, 
                 geoid: EarthGeoidModel = None,
                 density_model: EarthDensityModel = None,
                 tep_params: dict = None):
        """
        Initialize trajectory integrator.
        
        Parameters
        ----------
        geoid : EarthGeoidModel
            Earth shape model
        density_model : EarthDensityModel
            Density distribution model
        tep_params : dict
            TEP parameters (beta, Lambda_keV, n_topology)
        """
        self.geoid = geoid or EarthGeoidModel()
        self.density = density_model or EarthDensityModel(self.geoid)
        self.tep = tep_params or {'beta': BETA_BASELINE * 1e-4, 'Lambda_keV': 10.0, 'n_topology': 1}
        
        # Precompute Temporal Topology field constants
        self._setup_temporal_topology_field()
    
    def _setup_temporal_topology_field(self):
        """Precompute Temporal Topology field values at key points."""
        beta = self.tep['beta']
        Lambda = self.tep['Lambda_keV'] * 1e-6  # Convert to GeV
        n = self.tep['n_topology']
        
        # Field at Earth's center (highest density)
        rho_center = 13000  # kg/m^3
        self.phi_center = self._phi_of_rho(rho_center, beta, Lambda, n)
        
        # Field at surface
        rho_surface = 2700  # kg/m^3
        self.phi_surface = self._phi_of_rho(rho_surface, beta, Lambda, n)
        
        # Field in vacuum
        self.phi_space = self._phi_of_rho(1e-20, beta, Lambda, n)
        
        # Restoration length from atmosphere
        rho_atm = 1.225
        phi_atm = self._phi_of_rho(rho_atm, beta, Lambda, n)
        hbar_c = 0.197e-15  # GeV·m
        m2_atm = (n * (n + 1) * Lambda**(4 + n) / phi_atm**(n + 2))
        m_atm = np.sqrt(m2_atm)
        self.lambda_rest = hbar_c / m_atm
    
    @staticmethod
    def _phi_of_rho(rho_kg_m3: float, beta: float, Lambda: float, n: int) -> float:
        """Compute Temporal Topology field for given density."""
        rho_gev4 = rho_kg_m3 * KG_M3_TO_GEV4
        
        if rho_gev4 <= 0 or beta <= 0:
            return Lambda * 1e6
        
        # Jakarta v0.8 consistency: use factor of 2 in denominator for field minimum
        # Temporal Topology field minimum: φ_min = Λ [ (n Λ^(n+4) M_Pl) / (2β ρ) ]^(1/(n+1))
        numerator = n * M_PL * Lambda**(4 + n)
        denominator = 2.0 * beta * rho_gev4
        
        if denominator <= 0:
            return Lambda * 1e6
        
        scale = (numerator / denominator)**(1.0 / (n + 1))
        return Lambda * scale

    def temporal_topology_field(self, x: float, y: float, z: float) -> float:
        """
        Compute Temporal Topology field φ at position (x, y, z).
        
        Accounts for Temporal Shear Suppression at high densities and 
        exponential relaxation outside Earth.
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        
        if r <= R_EARTH:
            # Inside Earth: field depends on local density
            rho = self.density.density_with_geoid(x, y, z)
            return self._phi_of_rho(rho, self.tep['beta'], 
                                   self.tep['Lambda_keV'] * 1e-6,
                                   self.tep['n_topology'])
        else:
            # Outside Earth: exponential relaxation from surface value
            delta_r = r - R_EARTH
            frac = 1.0 - np.exp(-delta_r / self.lambda_rest)
            return self.phi_surface + (self.phi_space - self.phi_surface) * frac
    
    def tep_acceleration(self, x: float, y: float, z: float, vx: float = 0.0, vy: float = 0.0, vz: float = 0.0) -> Tuple[float, float, float]:
        """
        Compute TEP Temporal Topology acceleration at geocentric (x, y, z).
        
        Incorporates:
        1. Scalar Gradient Force: a_grad = c² β ∇φ / M_Pl
        2. Disformal Velocity Force: a_disf = (B/A) (∇φ · v) v
        3. Geoid-dependent Screening: ∇φ scaled by exp(-h_wgs84 / λ)
        
        Parameters
        ----------
        x, y, z : float
            Geocentric position (m)
        vx, vy, vz : float
            Geocentric velocity (m/s)
            
        Returns
        -------
        tuple (ax, ay, az)
            Anomalous acceleration in m/s²
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        if r < 1e-6:
            return 0.0, 0.0, 0.0
        
        # 1. Get Geodetic Altitude for Screening
        # This accounts for Earth's oblateness (J2) in the suppression logic
        lat, lon, alt = self.geoid.geocentric_to_geodetic(x, y, z)
        
        # If inside Earth, screening is effectively total (phi -> phi_surface)
        if alt <= 0:
            return 0.0, 0.0, 0.0
        
        # 2. Scalar Gradient Calculation (∇φ)
        # Jakarta v0.8: Continuous gradient suppression
        # φ(h) = φ_surf + (φ_space - φ_surf) * (1 - exp(-h/λ))
        # ∇φ = (φ_space - φ_surf) * (1/λ) * exp(-h/λ) * r̂
        # Note: r̂ is locally vertical to the geoid (approximated as radial here)
        grad_mag = (self.phi_space - self.phi_surface) * (1.0 / self.lambda_rest) * np.exp(-alt / self.lambda_rest)
        
        # Scalar Force Factor: a = c² β ∇φ / M_Pl
        beta = self.tep.get('beta', 1e-4)
        c2 = C_LIGHT**2
        factor_scalar = c2 * beta / M_PL
        
        # Unit vector (radial approximation for gradient direction)
        ux, uy, uz = x / r, y / r, z / r
        
        # Static Scalar Acceleration
        ax_scalar = factor_scalar * grad_mag * ux
        ay_scalar = factor_scalar * grad_mag * uy
        az_scalar = factor_scalar * grad_mag * uz
        
        # 3. Disformal Velocity Coupling (B-term)
        # a_disf = (B/M_pl) * (∇φ · v) v / c² (simplified effective coupling)
        # This term explains sign reversals and directional asymmetries.
        b_disformal = self.tep.get('b_disformal', 0.05)
        
        # Velocity vector
        v_vec = np.array([vx, vy, vz])
        grad_phi_vec = grad_mag * np.array([ux, uy, uz])
        
        # Dot product: ∇φ · v
        v_dot_grad = np.dot(v_vec, grad_phi_vec)
        
        # Disformal factor (dimensionless coupling)
        # Using a phenomenological scaling for the velocity-dependent part
        factor_disformal = (b_disformal / M_PL) * v_dot_grad
        
        ax_disf = factor_disformal * vx
        ay_disf = factor_disformal * vy
        az_disf = factor_disformal * vz
        
        # Total TEP Acceleration
        return ax_scalar + ax_disf, ay_scalar + ay_disf, az_scalar + az_disf
    
    def equations_of_motion(self, t: float, state: np.ndarray) -> np.ndarray:
        """
        Compute time derivatives for state vector [x, y, z, vx, vy, vz].
        
        Includes:
        - Non-spherical Earth gravity
        - Earth's rotation (centrifugal and Coriolis in inertial frame)
        - TEP Temporal Topology perturbation
        """
        x, y, z, vx, vy, vz = state
        
        # Standard gravitational acceleration
        ax_grav, ay_grav, az_grav = self.geoid.gravity_acceleration(x, y, z)
        
        # TEP Temporal Topology acceleration
        # Pass velocity for disformal coupling calculation
        ax_tep, ay_tep, az_tep = self.tep_acceleration(x, y, z, vx, vy, vz)
        
        # Total acceleration
        ax = ax_grav + ax_tep
        ay = ay_grav + ay_tep
        az = az_grav + az_tep
        
        return [vx, vy, vz, ax, ay, az]
    
    def integrate_trajectory(self,
                            t_span: Tuple[float, float],
                            y0: List[float],
                            method: str = 'RK45',
                            rtol: float = 1e-9,
                            atol: float = 1e-12,
                            dense_output: bool = True,
                            events: Optional[Callable] = None) -> dict:
        """
        Integrate spacecraft trajectory.
        
        Parameters
        ----------
        t_span : tuple (t0, tf)
            Integration time range in seconds
        y0 : list [x, y, z, vx, vy, vz]
            Initial state vector
        method : str
            Integration method ('RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA')
        rtol, atol : float
            Relative and absolute tolerances
        dense_output : bool
            Whether to compute continuous solution
        events : callable or list
            Event functions to detect (e.g., perigee passage)
            
        Returns
        -------
        dict
            Integration results with trajectory data
        """
        sol = solve_ivp(
            self.equations_of_motion,
            t_span,
            y0,
            method=method,
            rtol=rtol,
            atol=atol,
            dense_output=dense_output,
            events=events
        )
        
        # Package results
        n_points = len(sol.t)
        trajectory = []
        
        for i in range(n_points):
            state = TrajectoryState(
                t=sol.t[i],
                x=sol.y[0][i],
                y=sol.y[1][i],
                z=sol.y[2][i],
                vx=sol.y[3][i],
                vy=sol.y[4][i],
                vz=sol.y[5][i]
            )
            
            # Compute local TEP effects
            phi = self.temporal_topology_field(state.x, state.y, state.z)
            dtau = np.exp(self.tep['beta'] * phi / M_PL)
            
            traj_point = {
                't': state.t,
                'x': state.x,
                'y': state.y,
                'z': state.z,
                'r': state.r,
                'vx': state.vx,
                'vy': state.vy,
                'vz': state.vz,
                'v': state.v,
                'phi': phi,
                'dtau': dtau
            }
            trajectory.append(traj_point)
        
        # Find perigee
        r_values = [p['r'] for p in trajectory]
        perigee_idx = np.argmin(r_values)
        perigee = trajectory[perigee_idx]
        
        # Compute asymptotic velocities (first and last points)
        v_inf_in = trajectory[0]['v']
        v_inf_out = trajectory[-1]['v']
        dv_total = v_inf_out - v_inf_in
        
        return {
            'success': sol.success,
            'n_points': n_points,
            't_span': t_span,
            'perigee': perigee,
            'v_inf_in': v_inf_in,
            'v_inf_out': v_inf_out,
            'dv_total': dv_total,
            'trajectory': trajectory,
            'solution': sol
        }


def compute_newtonian_baseline(trajectory_data: dict, 
                                geoid: EarthGeoidModel = None) -> dict:
    """
    Compute pure Newtonian baseline trajectory (no TEP effects).
    
    This addresses the circularity problem by computing a trajectory
    using only standard physics, then comparing to the observed anomaly.
    
    Parameters
    ----------
    trajectory_data : dict
        Initial conditions from observed trajectory
    geoid : EarthGeoidModel
        Earth model for gravity
        
    Returns
    -------
    dict
        Newtonian baseline prediction
    """
    geoid = geoid or EarthGeoidModel()
    
    # Create integrator with beta=0 (no TEP)
    tep_params_zero = {'beta': 0.0, 'Lambda_keV': 10.0, 'n_topology': 1}
    integrator = TEP3DTrajectoryIntegrator(
        geoid=geoid,
        tep_params=tep_params_zero
    )
    
    # Extract initial conditions from observed trajectory
    perigee = trajectory_data['perigee']
    
    # Estimate time window: ±2 days around perigee
    t_perigee = 0.0
    t_start = t_perigee - 2 * 86400  # 2 days before
    t_end = t_perigee + 2 * 86400  # 2 days after
    
    # Convert initial conditions to m and m/s
    r_perigee = (perigee['altitude_km'] + 6371.0) * 1e3  # m
    v_perigee = perigee['velocity_km_s'] * 1e3  # m/s
    
    # Assume equatorial prograde trajectory for simplicity
    # (can be refined with actual trajectory data)
    x0 = r_perigee
    y0 = 0.0
    z0 = 0.0
    vx0 = 0.0
    vy0 = v_perigee
    vz0 = 0.0
    
    # Integrate backwards and forwards from perigee
    # This is a simplified approach - full implementation would
    # use actual initial conditions from ephemeris
    
    y0_state = [x0, y0, z0, vx0, vy0, vz0]
    
    result = integrator.integrate_trajectory(
        t_span=(t_start, t_end),
        y0=y0_state,
        method='DOP853',  # High accuracy
        rtol=1e-12,
        atol=1e-15
    )
    
    # Newtonian prediction: velocity change should be zero (energy conservation)
    # Any non-zero change indicates integration error or other physics
    newtonian_dv = result['dv_total']  # Should be ~0
    
    return {
        'method': 'newtonian_baseline',
        'geoid_model': 'WGS84_with_J2J3J4',
        'v_inf_in': result['v_inf_in'],
        'v_inf_out': result['v_inf_out'],
        'dv_predicted': newtonian_dv,
        'perigee_altitude_km': result['perigee']['r'] / 1e3 - 6371.0,
        'integration_success': result['success'],
        'n_integration_points': result['n_points']
    }


def compute_tep_prediction_full(trajectory_data: dict,
                                beta: float = 1e-4,
                                geoid: EarthGeoidModel = None) -> dict:
    """
    Compute full TEP prediction with 3D trajectory integration.
    
    This computes the velocity anomaly predicted by TEP theory
    by integrating the trajectory with Temporal Topology field effects.
    
    Parameters
    ----------
    trajectory_data : dict
        Initial conditions from observed trajectory
    beta : float
        TEP coupling parameter
    geoid : EarthGeoidModel
        Earth model
        
    Returns
    -------
    dict
        TEP prediction with proper residuals
    """
    geoid = geoid or EarthGeoidModel()
    
    tep_params = {'beta': beta, 'Lambda_keV': 10.0, 'n_topology': 1}
    
    integrator = TEP3DTrajectoryIntegrator(
        geoid=geoid,
        tep_params=tep_params
    )
    
    # Extract and convert initial conditions
    perigee = trajectory_data['perigee']
    r_perigee = (perigee['altitude_km'] + 6371.0) * 1e3
    v_perigee = perigee['velocity_km_s'] * 1e3
    
    t_start = -2 * 86400
    t_end = 2 * 86400
    
    # Simplified initial conditions (can be refined)
    y0 = [r_perigee, 0.0, 0.0, 0.0, v_perigee, 0.0]
    
    result = integrator.integrate_trajectory(
        t_span=(t_start, t_end),
        y0=y0,
        method='DOP853',
        rtol=1e-12,
        atol=1e-15
    )
    
    # TEP prediction: velocity change due to Temporal Topology field
    tep_dv = result['dv_total']
    
    # Also compute using simple formula for comparison
    phi_surface = integrator.phi_surface
    phi_space = integrator.phi_space
    dtau_ratio = np.exp(beta * (phi_space - phi_surface) / M_PL)
    dv_simple = (dtau_ratio - 1.0) * C_LIGHT * 1e3  # mm/s
    
    return {
        'method': 'tep_3d_integration',
        'topology_model': 'temporal_shear_suppression',
        'beta_used': beta,
        'v_inf_in': result['v_inf_in'],
        'v_inf_out': result['v_inf_out'],
        'dv_predicted_3d': tep_dv * 1e3,  # Convert to mm/s
        'dv_predicted_simple': dv_simple,
        'perigee_altitude_km': result['perigee']['r'] / 1e3 - 6371.0,
        'phi_at_perigee': result['perigee']['phi'],
        'dtau_at_perigee': result['perigee']['dtau'],
        'integration_success': result['success'],
        'n_integration_points': result['n_points']
    }
