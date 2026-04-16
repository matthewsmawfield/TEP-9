"""
Full 3D Field Integration for TEP Scalar Force

This module implements full 3D integration of the scalar field gradient
along flyby trajectories, replacing the simplified scalar force formula.

The scalar force is computed as:
F_φ = -β ∇φ/m

Where ∇φ is computed numerically from the full 3D field solution.

Key features:
- Numerical field gradient calculation
- Full trajectory integration
- 3D density profile integration
- Comparison with simplified formula
"""

import numpy as np
import json
from pathlib import Path
import sys
from scipy.integrate import quad
from scipy.interpolate import RegularGridInterpolator

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.step_logger import StepLogger


class Field3DIntegrator:
    """
    3D field integration for TEP scalar force calculation.
    
    Computes the velocity anomaly by integrating the scalar force
    along the full 3D trajectory.
    """
    
    # Physical constants
    M_PL = 2.435e18  # GeV - reduced Planck mass
    HBAR_C = 0.1973  # GeV·fm
    GEV_TO_J = 1.602e-10  # J/GeV
    C_LIGHT = 2.998e8  # m/s
    
    def __init__(self, beta=1e-4):
        """
        Initialize 3D field integrator.
        
        Args:
            beta: Coupling constant (dimensionless)
        """
        self.logger = StepLogger("step_017_3d_field_integration", PROJECT_ROOT)
        self.beta = beta
        
    def compute_field_gradient_3d(self, position, density_field):
        """
        Compute 3D gradient of scalar field at given position.
        
        Uses finite difference method on a 3D grid.
        
        Args:
            position: (x, y, z) tuple in meters
            density_field: Function returning density at position
            
        Returns:
            (∂φ/∂x, ∂φ/∂y, ∂φ/∂z) in GeV/m
        """
        delta = 1e3  # 1 km step size for finite difference
        
        x, y, z = position
        
        # Compute field values at neighboring points
        phi_x_plus = self.compute_field_at_position((x + delta, y, z), density_field)
        phi_x_minus = self.compute_field_at_position((x - delta, y, z), density_field)
        phi_y_plus = self.compute_field_at_position((x, y + delta, z), density_field)
        phi_y_minus = self.compute_field_at_position((x, y - delta, z), density_field)
        phi_z_plus = self.compute_field_at_position((x, y, z + delta), density_field)
        phi_z_minus = self.compute_field_at_position((x, y, z - delta), density_field)
        
        # Central difference
        dphi_dx = (phi_x_plus - phi_x_minus) / (2 * delta)
        dphi_dy = (phi_y_plus - phi_y_minus) / (2 * delta)
        dphi_dz = (phi_z_plus - phi_z_minus) / (2 * delta)
        
        return np.array([dphi_dx, dphi_dy, dphi_dz])
    
    def compute_field_at_position(self, position, density_field):
        """
        Compute scalar field value at given position using chameleon formula.
        
        φ = φ_min where φ_min = [(n+1)M_PlΛ^(4+n)/(ρβ)]^(1/(n+2))
        
        Args:
            position: (x, y, z) tuple in meters
            density_field: Function returning density at position
            
        Returns:
            Field value in GeV
        """
        rho = density_field(position)
        rho_gev4 = rho * 1.5e-41  # Convert kg/m³ to GeV^4
        
        # Chameleon parameters
        Lambda = 1e-3  # GeV
        n = 1
        
        if rho_gev4 <= 0 or self.beta <= 0:
            return Lambda * 1e6
        
        numerator = (n + 1) * self.M_PL * (Lambda ** (4 + n))
        denominator = rho_gev4 * self.beta
        
        if denominator <= 0:
            return Lambda * 1e6
        
        phi_min = (numerator / denominator) ** (1.0 / (n + 2))
        return phi_min
    
    def integrate_force_along_trajectory(self, trajectory_points, density_field, 
                                        velocity_profile):
        """
        Integrate scalar force along trajectory to compute velocity anomaly.
        
        Δv = ∫ F_φ/m dt = -β ∫ ∇φ/m dt
        
        Args:
            trajectory_points: List of (x, y, z) positions along trajectory
            density_field: Function returning density at position
            velocity_profile: Function returning velocity at each point
            
        Returns:
            Velocity anomaly in m/s
        """
        total_dv = np.array([0.0, 0.0, 0.0])  # dv_x, dv_y, dv_z
        
        for i in range(len(trajectory_points) - 1):
            pos1 = trajectory_points[i]
            pos2 = trajectory_points[i + 1]
            
            # Time step
            dt = 1.0  # seconds (simplified)
            
            # Compute field gradient at midpoint
            pos_mid = [(pos1[j] + pos2[j]) / 2 for j in range(3)]
            grad_phi = self.compute_field_gradient_3d(pos_mid, density_field)
            
            # Convert to SI units
            grad_phi_si = grad_phi * self.GEV_TO_J / 1e-15  # GeV/m to J/m
            
            # Scalar force per unit mass
            force_per_mass = -self.beta * grad_phi_si
            
            # Velocity increment
            dv = force_per_mass * dt
            total_dv += dv
        
        # Convert to mm/s
        dv_mm_s = total_dv * 1000
        
        return dv_mm_s
    
    def compute_simplified_force(self, trajectory_points, density_field):
        """
        Compute simplified scalar force (current implementation).
        
        Uses the formula: F = β ∇φ/m with ∇φ approximated from
        density gradient and screening.
        
        Args:
            trajectory_points: List of (x, y, z) positions
            density_field: Function returning density at position
            
        Returns:
            Velocity anomaly in mm/s
        """
        total_dv = 0.0
        
        for i in range(len(trajectory_points) - 1):
            pos1 = trajectory_points[i]
            pos2 = trajectory_points[i + 1]
            
            # Compute density gradient
            r1 = np.linalg.norm(pos1)
            r2 = np.linalg.norm(pos2)
            
            rho1 = density_field(pos1)
            rho2 = density_field(pos2)
            
            dr = r2 - r1
            dln_rho = np.log(rho2 / rho1) if rho2 > 0 and rho1 > 0 else 0
            
            # Simplified force
            if dr != 0:
                force = -self.beta * dln_rho / dr
                total_dv += force * dr
        
        return total_dv * 1000  # Convert to mm/s


class EarthDensityField:
    """
    3D density field for Earth including non-spherical components.
    """
    
    R_EARTH = 6.371e6  # m
    R_CORE = 3.48e6  # m
    J2 = 1.08263e-3  # Earth's oblateness parameter
    
    def __init__(self):
        self.logger = StepLogger("step_017_3d_field_integration", PROJECT_ROOT)
    
    def density_at_position(self, position):
        """
        Get density at 3D position including J2 oblateness.
        
        Args:
            position: (x, y, z) tuple in meters
            
        Returns:
            Density in kg/m³
        """
        x, y, z = position
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / r) if r > 0 else 0
        
        if r > self.R_EARTH:
            # Atmosphere
            altitude = r - self.R_EARTH
            scale_height = 8500
            rho = 1.225 * np.exp(-altitude / scale_height)
        elif r > self.R_CORE:
            # Mantle with J2 correction
            rho_base = 4500
            j2_correction = 1 + self.J2 * (3 * np.cos(theta)**2 - 1) / 2
            rho = rho_base * j2_correction
        else:
            # Core
            rho = 13000
        
        return max(rho, 1e-10)


def main():
    """Execute 3D field integration."""
    logger = StepLogger("step_017_3d_field_integration", PROJECT_ROOT)
    logger.section("STEP 017: 3D FIELD INTEGRATION")
    
    # Initialize integrator
    integrator = Field3DIntegrator(beta=1e-4)
    density_field = EarthDensityField()
    
    logger.subsection("INTEGRATION PARAMETERS")
    logger.info(f"β = {integrator.beta:.2e}")
    logger.info(f"J2 = {density_field.J2:.2e}")
    
    # Load trajectory data
    predictions_file = PROJECT_ROOT / 'results' / 'step004_tep_predictions.json'
    
    if not predictions_file.exists():
        logger.error("TEP predictions not found. Run step004 first.")
        logger.log_step_summary(0, "FAILED")
        return 1
    
    with open(predictions_file) as f:
        data = json.load(f)
    
    # Process each flyby
    logger.subsection("TRAJECTORY INTEGRATION")
    
    results = {}
    
    for name, pred in data['predictions'].items():
        if pred['observed']['dv_obs_mm_s'] != 0:
            logger.info(f"Processing {name}...")
            
            # Generate simplified trajectory (straight line through perigee)
            perigee_pos = np.array([0, 0, pred['perigee']['altitude_km'] * 1000])
            
            # Create trajectory points (simplified: radial approach)
            n_points = 100
            r_min = pred['perigee']['altitude_km'] * 1000
            r_max = r_min + 1e6  # 1000 km above perigee
            
            trajectory_points = []
            for r in np.linspace(r_max, r_min, n_points):
                trajectory_points.append((0, 0, r))
            for r in np.linspace(r_min, r_max, n_points):
                trajectory_points.append((0, 0, r))
            
            # Compute 3D integrated force
            dv_3d = integrator.integrate_force_along_trajectory(
                trajectory_points, density_field.density_at_position, None
            )
            
            # Compute simplified force
            dv_simplified = integrator.compute_simplified_force(
                trajectory_points, density_field.density_at_position
            )
            
            # Get current prediction
            dv_current = pred['tep_predictions']['dv_tep_mm_s']
            
            results[name] = {
                'dv_3d_mm_s': float(np.linalg.norm(dv_3d)),
                'dv_simplified_mm_s': float(dv_simplified),
                'dv_current_mm_s': dv_current,
                'ratio_3d_simplified': float(np.linalg.norm(dv_3d) / dv_simplified) if dv_simplified > 0 else 1.0
            }
            
            logger.info(f"  3D integration: {results[name]['dv_3d_mm_s']:.4f} mm/s")
            logger.info(f"  Simplified: {results[name]['dv_simplified_mm_s']:.4f} mm/s")
            logger.info(f"  Current: {results[name]['dv_current_mm_s']:.4f} mm/s")
            logger.info(f"  Ratio: {results[name]['ratio_3d_simplified']:.2f}")
    
    # Save results
    results_dir = PROJECT_ROOT / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = results_dir / 'step017_3d_field_integration_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.section("SAVING RESULTS")
    logger.info(f"Results saved to: {output_file}")
    
    logger.log_step_summary(0, "SUCCESS")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
