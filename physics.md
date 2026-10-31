# GeoFluxSim: physics and numerics

This page documents the model behind GeoFluxSim, version 1.0.0. The web app (`index.html`) and the Jupyter version (`notebooks/geofluxsim_jupyter.py`) implement the same equations.

## 1. Setting

- Flume length L = 8 m, width W = 0.6 m, discretised into 240 cells (Δx ≈ 3.3 cm).
- The bed descends downstream at slope S: z_b(x) = S (L − x).
- A current of thickness h₀, velocity U₀ and volume concentration C₀ enters at x = 0 for a set duration. Afterwards the inlet is a closed wall. The outlet is open.

## 2. Sediment

Six grain classes: siliciclastic and carbonate, each with fine, median and coarse sizes

D = D₅₀ σ_g⁻¹, D₅₀, D₅₀ σ_g

carrying 25%, 50% and 25% of that component's sediment. σ_g is the sorting. The carbonate fraction splits the inflow concentration between the two components.

Submerged specific gravity of class i: Rᵢ = (ρ_s,i − ρ_a) / ρ_a, where ρ_a is the ambient water density.

## 3. Governing equations (layer-averaged; Parker et al., 1986)

Fluid volume:

∂h/∂t + ∂(Uh)/∂x = e_w |U|

Momentum:

∂(Uh)/∂t + ∂(U²h + ½ g′h²)/∂x = g′ h S − C_d U|U|,  with g′ = g Σ Rᵢ Cᵢ

Sediment, for each class i:

∂(Cᵢh)/∂t + ∂(U Cᵢ h)/∂x = − r₀ w_s,i (1 − C)^4.65 Cᵢ

### Water entrainment (Parker et al., 1987)

e_w = a / √(1 + 718 Ri^2.4),  Ri = g′h / U²,  a = 0.075 by default

### Settling velocity (Ferguson & Church, 2004)

w_s = R g D² / (C₁ ν + √(0.75 C₂ R g D³))

Shape constants (C₁, C₂):

| Grain shape | C₁ | C₂ |
|---|---|---|
| Rounded quartz | 18 | 1.0 |
| Angular quartz | 20 | 1.2 |
| Platy, mica-rich | 24 | 1.5 |
| Ooids / peloids (near-spherical) | 18 | 0.6 |
| Bioclastic fragments (angular) | 20 | 1.2 |
| Plates & rods (skeletal) | 24 | 1.5 |

Kinematic viscosity from water temperature T (°C): ν = 1.79×10⁻⁶ / (1 + 0.03368 T + 0.000221 T²) m²/s.

Hindered settling uses the Richardson & Zaki (1954) factor (1 − C)^4.65, where C is the total volume concentration.

### Deposition

The deposition rate of class i is r₀ w_s,i (1 − C)^4.65 Cᵢ, where r₀ is the ratio of near-bed to layer-averaged concentration (default 1.6). Deposited sediment becomes bed thickness through the porosity λ (default 0.4):

thickness = deposited volume / (1 − λ)

Deposits are stored in 1.5 s layers, which gives the stratigraphy, time lines, cores and downstream distributions.

## 4. Vertical structure (for profiles and the 3D view)

The layer model is depth-averaged. Profiles at the probe are reconstructed as follows.

**Velocity** (after Altinakar et al., 1996), with η = z/h:
- below the maximum at η_m = 0.2: u ∝ (η/η_m)^(1/6)
- above it: u ∝ exp(−1.5 ((η − η_m)/(1 − η_m))²)

The profile is scaled so that ∫u dz = Uh.

**Concentration:** c(η) = C · k/(1 − e^(−k)) · e^(−kη), with a smooth cut-off above the layer top. k is chosen so that the near-bed value equals r₀ C.

## 5. Numerical method

- Finite-volume scheme with a Rusanov (local Lax–Friedrichs) flux, wave speed |U| + √(g′h).
- Explicit time step limited by CFL = 0.4 (and at most 0.05 s).
- Velocity is desingularised where the layer is vanishingly thin (h < 1 mm). The minimum layer thickness is 10⁻⁴ m.
- Bed friction is treated implicitly. Deposition uses an exact exponential decay over each step, so concentrations stay non-negative.
- The front is the furthest cell where Cᵢh exceeds 0.1% of the inflow value.
- Sediment budget: injected = deposited + bypassed through the outlet + still in suspension. It closes to within rounding.

## 6. Assumptions and limitations

- Width-averaged: no lateral spreading. The deposit is uniform across the flume.
- No bed erosion or resuspension.
- Single layer: no separate head dynamics, internal stratification or hydraulic-jump modelling beyond what the shallow-water equations capture.
- The 3D turbulent billows are procedural, for visual realism. They follow the computed thickness and concentration but are not simulated turbulence.
- Default coefficients are typical literature values. They have not been calibrated against specific experiments.

## References

- Altinakar, M. S., Graf, W. H., & Hopfinger, E. J. (1996). Flow structure in turbidity currents. *Journal of Hydraulic Research*, 34(5), 713–718.
- Ferguson, R. I., & Church, M. (2004). A simple universal equation for grain settling velocity. *Journal of Sedimentary Research*, 74(6), 933–937.
- Parker, G., Fukushima, Y., & Pantin, H. M. (1986). Self-accelerating turbidity currents. *Journal of Fluid Mechanics*, 171, 145–181.
- Parker, G., Garcia, M., Fukushima, Y., & Yu, W. (1987). Experiments on turbidity currents over an erodible bed. *Journal of Hydraulic Research*, 25(1), 123–147.
- Richardson, J. F., & Zaki, W. N. (1954). Sedimentation and fluidisation: Part I. *Transactions of the Institution of Chemical Engineers*, 32, 35–53.
