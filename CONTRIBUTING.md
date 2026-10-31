# Contributing to GeoFluxSim

Thank you for helping improve GeoFluxSim. Feedback from people using it in class, in the lab or in industry is the most valuable contribution.

## Reporting a problem or suggesting an idea
[Open an issue](../../issues/new/choose) and pick the template that fits:
- **Bug report:** something doesn't work or looks wrong.
- **Feature idea:** something you'd like it to do.
- **Classroom or research feedback:** how you used it and what would make it more useful.

Please include your browser and device, the settings you used, and a screenshot if you can.

## Suggesting a code change
1. Fork the repository and create a branch.
2. Make your change in `index.html`. If you change the physics, mirror it in `notebooks/geofluxsim_jupyter.py` and update `docs/physics.md`.
3. Run through `docs/release-checklist.md`.
4. Open a pull request that explains what changed and why.

## Changes to the physics
Please cite the source of any new equation or coefficient. Where possible, show that the result still closes the sediment budget and behaves sensibly against the sanity checks in the release checklist.
