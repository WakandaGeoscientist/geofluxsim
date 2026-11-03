2026-11 update

Added

	•	Study builder: vary one or two settings, compare setups, or test which settings matter most; choose the results to show, including vertical profiles and grain-size distributions; export as CSV.
	•	Outcrop and synthetic seismic views with AVO analysis (angle gathers, AVO curves, intercept–gradient crossplot, brine/oil/gas).
	•	"How the flow built the deposit": bypass, erosion and deposition zones against the Froude number.
	•	Mud types (kaolinite, illite, chlorite, bentonite, carbonate mud, silt-like), several per run, with cohesion after Baker et al. (2017, 2025).
	•	Deposit photo reader and hypothesis testing; comparison videos; model card.

Changed

	•	Default experiment: a steep upper slope easing to a flat basin floor, so the current bypasses the slope and builds a lobe.
	•	Tidier interface: foldable sections, a shortcut row on each tab, a slimmer run bar, Simple and Advanced levels.
	•	Share links point to whichever website hosts the app.

Fixed

	•	Compare runs: the "differs by" column for pinned runs.


# Changelog

All notable changes to GeoFluxSim. Versions follow MAJOR.MINOR.PATCH:
- MAJOR for changes to the model that alter results;
- MINOR for new features;
- PATCH for fixes.

## [1.0.0] – 2026-10-31
First public release.

### Added
- Layer-averaged turbidity-current model with six grain classes (siliciclastic and carbonate × fine, median, coarse).
- Ferguson & Church settling with grain shape options, temperature-dependent viscosity and hindered settling.
- 3D flume with realistic, concentration and velocity views, cut-away centre slice and probe profile.
- Layered deposit with time lines, cross-section, core, downstream distribution and summary statistics.
- Live physics diagram and equation list at the probe, settling-velocity chart and sediment budget.
- Adjustable sediment, fluid, flow and model coefficients.
- Jupyter version with the same physics and plots.
