# WCSim pion-analysis output — conventions

Conventions for the custom trees added by this WCSim modification
(`pion_photons`, `secondary_tracks`, `pion_steps`) and the things you must know
to read them correctly. Ship this file with the modified source.

## Units and frame

- Positions: **cm**. Times: **ns**. Energies: **MeV** (kinetic energy unless a
  branch name says otherwise).
- Coordinate frame: the WCSim detector frame (same frame as the stock `track_*`
  truth tree). Positions are `GetPosition()` in that frame, converted to cm.
- Optical-photon energies (`pp_*` start/end energy, if present) are ~3e-6 MeV.
  They are **filled**, not zero — they only round to 0 in a quick printout.
  Convert to eV before displaying.

## `src_cat` legend

Used by `pp_src_cat`, `st_src_cat`, `pst_src_cat` (the integer is stored so you
can cut on it directly):

| value | name          | meaning                                             |
|------:|---------------|-----------------------------------------------------|
| 0 | PrimaryPion   | primary π± (Geant4 parentID == 0)                       |
| 1 | SecPionPlus   | secondary π+ (e.g. from hadronic inelastic)             |
| 2 | SecPionMinus  | secondary π−                                            |
| 3 | PionZero      | π⁰, or a particle in the π⁰→γ→e± chain                   |
| 4 | Muon          | μ± (from π decay or any source)                         |
| 5 | MichelElec    | e± from μ decay ("Decay", parent is μ±)                 |
| 6 | DeltaRay      | e± from ionisation (hIoni/eIoni/muIoni/…)               |
| 7 | Nuclear       | nuclear-breakup fragment: p, n, α, heavy recoil         |
| 8 | Gamma         | γ (nCapture, bremsstrahlung, unmatched π⁰, …)           |
| 9 | Other         | anything else                                           |

(Definition: `SourceCategory` enum in `include/WCSimAncestryMap.hh`.)

## The custom trees

npz prefixes below follow the analysis-side naming (`pp_`, `st_`, `pst_`); the
ROOT tree names are `pion_photons`, `secondary_tracks`, `pion_steps`.

### `pion_photons` (pp_) — one row per Cherenkov photon that made a PE
- `pp_edir_{x,y,z}` is the photon emission direction, a proper **unit vector**
  (|edir| = 1). The PMT incidence angle is **not** stored — recompute it offline
  from `pp_edir`, `pp_hit_pmt`, and the PMT orientation in the geofile.
- `pp_creator` is always `Cerenkov` and `pp_end_process` always `Transportation`
  by construction — these are correctly constant, not bugs.

### `secondary_tracks` (st_) — one row per interesting secondary track
Interesting = π±, π⁰, μ±, Michel/conv e±, π⁰→γ, and nuclear fragments
(`SecTrack_IsInteresting` in `src/WCSimSecondaryTrackTree.cc`).
- **One row per track.** `x0/y0/z0/ke0` = birth (track vertex);
  `x1/y1/z1/ke1` = end. `(x1−x0, …)` is start-to-end displacement, not a path
  length.
  Geant4 suspends and resumes a track many times, calling
  Pre/PostUserTrackingAction on each segment. Earlier output therefore held
  one row per *segment* (486 contiguous rows for a single 235 MeV muon, with
  `x0` pinned at the track origin while `ke0` tracked the segment start).
  Birth is now cached only on first sight and a row is emitted only once the
  track reaches a terminal status, so the row count is the track count.
- `st_parent_trk` is the **real Geant4 parent track id** (use this for
  parentage). `st_parent_pdg` is the parent PDG.
- `st_end_process` is the true terminating process, including at-rest death
  (see the A1/A2 fix note below).

### `pion_steps` (pst_) — one row per pion/muon interaction or decay/stop vertex
Records only discrete vertices of π/μ tracks, not every ionisation step:
- in flight: `hadElastic`, `pi{+,-,0}Inelastic`, `Decay`, and hard `CoulombScat`
  (muons);
- **at rest** (after the A1 fix): the vertex where a π/μ ranges out and dies.
  These have `ke_pre ≈ 0`. Rows with process `StopAndKill_AtRest` mean the
  track stopped but no terminating process could be resolved (should be rare —
  investigate if common).

  Process names actually produced by the current physics list, measured on
  100-event runs and confirmed on a 10,000-event pi+ @198 MeV production:

  | process | particle | meaning |
  |---|---|---|
  | `Decay` | π+ | π+ → μ+ν at rest |
  | `Decay` | μ+ | Michel decay at rest |
  | `hBertiniCaptureAtRest` | π− | π− nuclear capture at rest |
  | `muMinusCaptureAtRest` | μ− | μ− capture **and** decay-in-orbit (Geant4 reports both under this one name) |

  NOTE: π− capture is named **`hBertiniCaptureAtRest`**, not `pi-CaptureAtRest`.
  The `isInteresting` list in `WCSimSteppingAction.cc` still names
  `pi±CaptureAtRest`, which this physics list never emits; at-rest rows are
  recorded anyway because at-rest death is captured regardless of process name.
  Cut on `ke_pre ≈ 0` rather than on a process-name list.
- `n_secondaries` = **physics daughters** produced at that vertex. Optical
  (Cherenkov/scintillation) photons are excluded; counting them made this
  branch scale with the emitter's velocity rather than count daughters
  (e.g. 73 at a CoulombScat vertex, which has no daughters at all).

#### Null inelastic vertices — filter these out

A `pi+Inelastic` row with `n_secondaries == 0` is **not a real interaction** and
must be excluded when counting inelastic scatters. Measured on 10,000 pi+ @198
MeV events:

| vertex | continues the track | ends the track |
|---|---:|---:|
| `pi+Inelastic`, `n_secondaries > 0`  |    0 | 6223 |
| `pi+Inelastic`, `n_secondaries == 0` |   50 |    0 |

The separation is exact. A real inelastic interaction always destroys the pion
and spawns daughters; these always leave the pion alive with no daughters. They
are null/failed interactions — Geant4 picked `pi+Inelastic` as the step-limiting
process but the final-state generator produced nothing and returned the
projectile intact. They are 0.8% of inelastic vertices.

```python
# when counting inelastic interactions
real_inelastic = (process == "pi+Inelastic") & (n_secondaries > 0)
```

By contrast, **`hadElastic` with `n_secondaries == 0` is normal and physical**
(715 of 766 zero-daughter vertices in that run): the nuclear recoil falls below
Geant4's production threshold, so no secondary track is created while the pion
simply changes direction. Keep those rows.

This was only diagnosable once optical photons were removed from
`n_secondaries`; under the old padded counting these vertices were invisible.

## The stock `track_*` truth tree — NOT part of this modification

`track_*` is WCSim's original `jhfNtuple` geant tree. Its quirks are upstream
WCSim behaviour (see `WCSimEventAction.cc` ~L895–970), not something this
modification introduces — do not "fix" them in the sim; handle them in analysis:

- **Generator ghost.** `track_flag == -1` is the "incoming neutrino" row, given
  an imaginary start ~100 m upstream (`start_z ≈ -10000 cm`). Remove it with
  `start_z > -200` (cm).
- **Null/target row.** `track_flag == -2` is the "target" row (pdg 0, zero
  fields). Filter it out.
- Because of those two rows, `track_id == 0` appears twice per event.
- **`track_parent` is a PDG code, not a track id** (e.g. 211, −13, 999 sentinel).
  For real parentage use `secondary_tracks.st_parent_trk` instead.

## Known open items (as of this writing)

- `digi_hit_trigger`, `trigger_time` are all-zero and `trigger_type` constant:
  the current productions are effectively untriggered (DAQ/trigger config), not a
  tree bug. Turn triggering on in the `.mac` if you need it.
- π⁻: the capture-at-rest path now works and is exercised — a 10k pi+ run
  yields 28 `hBertiniCaptureAtRest`, 5 `pi-Inelastic` and 2 `hadElastic` π⁻
  vertices from secondary pions. A dedicated stopping-π⁻ production is still
  needed for real statistics; earlier π⁻ files predate the modification and
  have no `pion_steps` tree at all.
- `G4Navigator` "Track stuck or not moving" warnings appear at ~0.3% of events
  in volume `-CDS-`. This is pre-existing WCTE geometry behaviour, unrelated to
  these trees.

## Production check

Run `check_constant_branches.py` once per production and fail loudly:

```
python3 pion_analysis_tools/check_constant_branches.py /path/to/output_dir
```

It flags any branch that is silently constant or empty, with an allow-list for
the ones that are legitimately constant. For **mixed-species or multi-energy**
runs, remove `pid` and `energy` from `ALLOW_CONSTANT` so they are required to
vary.
