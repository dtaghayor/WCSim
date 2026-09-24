# WCSim pion-analysis modification

A fork of [WCSim](https://github.com/WCSim/WCSim) adding truth-level bookkeeping
for charged-pion interactions, aimed at π⁺/π⁻ topology and Michel-electron
studies in WCTE. Branched from upstream `develop`.

Three ROOT trees are added alongside the standard WCSim output, plus a full
track-ancestry map so every Cherenkov photon can be traced back to the physics
process that ultimately produced it.

| tree | one row per | purpose |
|---|---|---|
| `pion_steps` | π/μ interaction, decay or stop vertex | where pions scatter, interact, decay or come to rest |
| `secondary_tracks` | interesting secondary track | π±, π⁰, μ±, Michel/conv e±, nuclear fragments, with birth and death kinematics |
| `pion_photons` | Cherenkov photon that produced a PE | emission point/direction/wavelength plus the ancestry category of its source |

Every row carries a `src_cat` code identifying the physics origin of the
particle (primary pion, secondary π±, π⁰ chain, muon, Michel electron, delta
ray, nuclear fragment, gamma, other).

**Read [`OUTPUT_CONVENTIONS.md`](OUTPUT_CONVENTIONS.md) before using the output.**
It documents the units and frame, the `src_cat` legend, the exact semantics of
each tree, which vertices must be filtered out, and which quirks belong to
upstream WCSim rather than to this modification.

## What changed relative to upstream

New files:

```
include/WCSimAncestryMap.hh        src/WCSimAncestryMap.cc
include/WCSimPionStepsTree.hh      src/WCSimPionStepsTree.cc
include/WCSimSecondaryTrackTree.hh src/WCSimSecondaryTrackTree.cc
include/WCSimPionSourceTree.hh     src/WCSimPionSourceTree.cc
pion_analysis_tools/
```

Modified: `WCSimSteppingAction`, `WCSimTrackingAction`, `WCSimEventAction`,
`WCSimRunAction`, `WCSimWCSD`, and `CMakeLists.txt`.

Sources are picked up by the existing `file(GLOB ...)` in `src/CMakeLists.txt`,
so no build-system changes are needed to compile the additions.

### Correctness fixes worth knowing about

These were found by auditing the output against the code and are the reason the
trees can be trusted:

- **At-rest termination is recorded.** Geant4 reports the *along-step* process
  (often `Scintillation`) for a track that dies at rest, not the process that
  killed it. A stopping pion's decay and a π⁻/μ⁻ capture were therefore missing
  entirely. They are now detected via `fStopAndKill` with zero pre-step kinetic
  energy, and the true process is recovered from the creator process of the
  non-optical secondaries. On a 10k π⁺ run this recovers 13,262 at-rest vertices
  that previously did not appear at all, and takes the fraction of events with
  at least one recorded vertex from ~90% to 100%.
- **`secondary_tracks` is one row per track.** Geant4 suspends and resumes
  tracks, calling `Pre`/`PostUserTrackingAction` on every segment, which
  previously produced one row per segment (486 rows for a single 235 MeV muon).
  Birth is now cached once and a row emitted only at a terminal status.
- **`n_secondaries` counts physics daughters**, excluding Cherenkov and
  scintillation photons. Previously it scaled with the emitter's velocity rather
  than counting daughters.

## Building

Standard WCSim build; it needs ROOT and Geant4 (developed against Geant4 10.3.3
and ROOT 6.28):

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=install
make install -j4
source install/bin/this_wcsim.sh     # sets WCSIM_BUILD_DIR, required at runtime
```

`WCSIM_BUILD_DIR` must be set or WCSim exits at startup.

## Running

```bash
WCSim <steering_macro> <tuning_parameter_macro>
```

The added trees are written unconditionally; no macro commands are required to
enable them.

## Validation

Checked against a 10,000-event π⁺ at 198 MeV in the WCTE geometry, plus
100-event μ⁺, μ⁻ and π⁺ runs. In-flight interaction rates agree with the
unmodified build (`hadElastic` 71.8 vs 71.5 per 100 events; `pi+Inelastic` 62.7
vs 61.8), confirming the additions do not perturb the physics. `pion_steps`
recovers μ⁻ capture-at-rest and π⁻ capture-at-rest, neither of which appeared
in earlier output.

`pion_analysis_tools/check_constant_branches.py` is a per-production sanity
check that fails loudly on any branch which is silently constant or empty — the
failure mode that hid several unfilled branches for a long time:

```bash
python3 pion_analysis_tools/check_constant_branches.py /path/to/output_dir
```
