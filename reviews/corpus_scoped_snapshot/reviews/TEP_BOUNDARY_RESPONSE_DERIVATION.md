# TEP boundary response: first constructive derivation

19 September 2026. Objective: derive a physical prediction from the canonical matter metric that can eventually be tested against observations. This note establishes a conditional response theorem within the companion programme, whose field-existence premise is imported for TEP-9. The new inference is the particular outer-system configuration and its observational response.

## 1. Start from the canonical action for matter

The repository's `core/definitions.md` specifies the matter metric

\[
\widetilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}
 +B(\phi)\partial_\mu\phi\partial_\nu\phi.
\]

Take the conformal limit B=0 and write psi=ln A. In a local patch where the background metric can be approximated by Minkowski space, the action of an unscreened structureless test particle is

\[
S_p=-mc^2\int e^{\psi(\boldsymbol x)}
       \sqrt{1-v^2/c^2}\,dt.
\]

This is a derived consequence of minimal coupling to the stated matter metric. It does not add an independent force law. The patch approximation isolates the boundary contribution; a Solar-System application must restore the background gravitational field, any disformal terms, and a justified response for extended screened bodies.

Consistency of equations under conformal frame changes is a standard requirement; see [Morris, Physical Review D 90, 107501 (2014)](https://doi.org/10.1103/PhysRevD.90.107501). The boundary formulas below follow directly from the displayed particle action, rather than from a fit to the comet data.

## 2. Static planar boundary: an exact conditional result

Let psi depend only on the normal coordinate z. The action is invariant under translations of t and both coordinates parallel to the boundary. With gamma=(1-v²/c²)^(-1/2), the conserved energy and parallel momentum are

\[
E=mc^2 A\gamma,\qquad
\boldsymbol p_\parallel=mA\gamma\boldsymbol v_\parallel.
\]

Consequently the parallel coordinate velocity is constant. For a transmitted particle crossing between constant asymptotic factors A_in and A_out,

\[
v_{n,\mathrm{out}}^2
=c^2-v_\parallel^2
 -\left(\frac{A_\mathrm{out}}{A_\mathrm{in}}\right)^2
  (c^2-v_\mathrm{in}^2).
\]

The positive square root describes transmission in the original direction. If the required square becomes negative, this transmitted branch is inaccessible; turning/reflection must be treated instead.

For small delta psi and nonrelativistic motion this becomes

\[
v_{n,\mathrm{out}}^2-v_{n,\mathrm{in}}^2
 \simeq -2c^2\Delta\psi.
\]

For a perturbatively small change relative to the initial normal kinetic energy,

\[
\Delta v_n\simeq-\frac{c^2\Delta\psi}{v_{n,\mathrm{in}}}.
\]

These are coordinate-velocity relations in the specified background patch. Observed astrometry requires the local matter-frame clocks, observer motion and photon measurement model. The velocity angle change is not itself a heliocentric eccentricity-vector or periapsis angle.

## 3. Clock response and trajectory response have the same amplitude

Expanding the particle action gives, up to the leading retained order,

\[
L=\tfrac12mv^2-m\Phi_N-mc^2\psi+\mathrm{constant},
\qquad
\boldsymbol a=-\boldsymbol\nabla\Phi_N-c^2\boldsymbol\nabla\psi.
\]

On a specified worldline the corresponding conformal clock factor is d(tau_tilde)=A d(tau_g). Thus the fractional factor between asymptotic regions is exp(delta psi)-1. A model cannot independently choose the lapse change and the trajectory impulse in this limit. This connection supplies a route to a joint orbital/clock consistency test.

It also distinguishes two geometries that a generic phrase such as "domain boundary" leaves ambiguous:

- A step with different asymptotic A values can change the transmitted direction and normal kinetic energy.
- A static, planar, localized bump with equal asymptotic A values returns a transmitting particle to its original asymptotic velocity. The transit time and position along the boundary can still change. Curved, time-dependent or disformal fields need a separate calculation; this cancellation is not a no-go theorem for them or for all TEP models.

The result is constructive: TEP's conformal matter sector can produce a calculable directional response. Its sign, energy change and clock response are linked, so the explanation is testable rather than an unrestricted re-description of an angle.

## 4. Numerical check of the response, not fabricated evidence

`scripts/experimental/conformal_boundary_response.py` integrates the leading equation through a tanh step and a sech² bump at ten signed/zero amplitudes. These amplitudes and initial velocities are deliberately prescribed synthetic controls, not estimated Solar-System parameters. It checks the analytic energy relation, zero-field limit, equal-endpoint cancellation and tighter-tolerance convergence. Results are in `results/experiments/conformal_boundary_response.json`.

The tanh and sech² profiles are test functions. The repository's local `core/scalar_field.py` implements a laboratory logarithmic ansatz; it does not supply a solved heliospheric boundary field. Nothing in this calculation promotes a laboratory fit into an outer-Solar-System prediction or establishes boundary stability or scalar stress-energy backreaction.

## 5. The empirical link to complete

Use a physical field solution, with its parameters fixed or estimated on explicitly designated training observations, to compute trajectories and synthetic measurements. Refit those measurements with the same conventional nuisance model used for the real astrometry. Then compare the surviving predicted signal to the real observation residuals, with covariance and survey selection accounted for. The same field must satisfy its clock and energy predictions.

The most useful next mathematical object is the derivative of those fitted measurement residuals with respect to the field parameters. It turns this conditional boundary response into an identifiable signal model. The field profile, its scalar dynamics and screening law must be specified before that derivative can be claimed as the canonical TEP prediction. Inferring them from the desired orbital answer would not complete the proof.
