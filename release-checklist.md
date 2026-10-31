# Release checklist

Run through this before pushing a new version to `main`. The push publishes the live site.

## Behaviour
- [ ] Page loads with no errors in the browser console (F12 → Console).
- [ ] Play, Pause and Reset work. The run ends cleanly at 90 s and shows "Run complete".
- [ ] Changing any parameter restarts the run. The probe, display and speed settings do not.
- [ ] Camera: 3D overview, Side wall and Probe camera all render. Drag to orbit and zoom work.
- [ ] Flow views: Realistic, Concentration and Velocity all render, with the colour bar and white probe profile.
- [ ] Layer toggles: Flow and Deposit can each be hidden.
- [ ] All plots update while playing: profiles, front, budget, deposit cross-section, core, distribution, physics diagram, settling chart and equation list.
- [ ] Deposit colours switch between Composition and Grain size.

## Physics sanity checks (default settings, run to the end)
- [ ] Deposited ≈ 75% and bypassed ≈ 25% of injected sediment.
- [ ] Deposit is thickest at the inlet (about 7 mm) and thins downstream.
- [ ] Coarser grains give a thicker, shorter deposit. Finer grains carry further.

## Devices
- [ ] Chrome or Edge on desktop
- [ ] Safari on iPhone or iPad
- [ ] Chrome on Android
- [ ] Light and dark mode

## Housekeeping
- [ ] Version number updated in `CHANGELOG.md` and `CITATION.cff`.
- [ ] If physics changed: `docs/physics.md` updated and the Jupyter version kept in step.
- [ ] After publishing, check https://geofluxsim.netlify.app in a private window.
