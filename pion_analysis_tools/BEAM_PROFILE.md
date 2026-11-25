# WCTE beam profile for the GPS generator

These commands make GPS beam events look like the real WCTE beam: a spot that
has passed through T5, and a momentum spread set by the collimators. Each one
is **off by default**, so existing macros behave exactly as before.

The spot and the momentum spread come from Loris's and Thorsten's work in
[Martiloris/WCSim](https://github.com/Martiloris/WCSim). The T5-shaped
acceptance and the momentum-space spread were added here.

## Commands

Put these after the `/gps/...` lines:

```
/mygen/gpsBeamT5Acceptance true      # keep only particles that passed through a T5 bar
/mygen/gpsBeamSigmaX 5 cm            # Gaussian beam width across the beam (x)
/mygen/gpsBeamSigmaY 5 cm            # ... and (y)
/mygen/gpsBeamMomentumSpread 0.02    # sigma_p / p (0.01 or 0.02, from the collimator setting)
```

| Command | Default | What it does |
|---|---|---|
| `/mygen/gpsBeamSigmaX`, `/mygen/gpsBeamSigmaY` | 0 | Width of a 2D Gaussian spread of the start point across the beam, centred on `/gps/position`. |
| `/mygen/gpsBeamT5Acceptance` | false | Rejects and redraws any start point whose offset misses all 8 T5 bars. |
| `/mygen/gpsBeamPipeDiameter` | 0 (off) | Optional circular cut instead of, or as well as, T5. This is the original Martiloris option. |
| `/mygen/gpsBeamMomentumSpread` | 0 (off) | Redraws each primary's \|p\| from Gauss(p0, f·p0) and keeps its direction. p0 comes from `/gps/energy`, which is kinetic energy. |

The run log states which options are active once per run, for example:

```
GPS momentum spread ON: sigma_p/p = 0.02, p0 = 307.366 MeV/c (T0 = 198 MeV), sigma_p = 6.14732 MeV/c, equivalent sigma_T ~ 5.59728 MeV
GPS beam spot ON: sigmaX 5 cm, sigmaY 5 cm, acceptance: T5-bars, axes u=(1,0,0) v=(0,1,0)
```

### How the spot works

- The spread is applied in the plane perpendicular to the beam direction. For a
  beam along ±z, σx and σy act along global x and y.
- T5 is centred on the `/gps/position` axis. Each of the 8 bars is 16.25 mm high,
  stacked along y from bar 1 (bottom) to bar 8, with lengths along x of
  41/94/112/123/123/112/94/41 mm. The bar table is in
  `include/WCSimT5Geometry.hh`, which the geometry also uses, so the acceptance
  and the T5 volumes always match.
- Only the start point's x and y change. z stays at `/gps/position`. The
  spread angle of the beam between T5 and the start point is ignored.
- The momentum spread is Gaussian in **momentum**, not kinetic energy, so you
  no longer need to convert σ_T = (p²/E)·σ_p/p by hand for `/gps/ene/type Gauss`.

## Why σ = 5 cm with the T5 shape

WCTE data on T5 hit positions (2714 muons; pions are consistent) show a beam
that is wide compared with T5. The best measure is the fraction of hits per
bar. It doesn't depend on how a bar number is turned into a position in mm.
Fitting a Gaussian beam with the T5 acceptance to those fractions gives
σ ≈ 4.5–5 cm, with bar 6 excluded:

| bar | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| simulation, σ = 5 cm | .039 | .114 | .154 | .194 | .184 | .151 | .121 | .043 |
| data (muons) | .036 | .115 | .178 | .196 | .223 | **.056** | .152 | .044 |

Bars 1, 2, 4 and 8 agree within the simulation's statistical error (about ±0.005).
Bar 6 has far fewer hits in data than its mirror, bar 3. A symmetric beam can't
produce that, so it is most likely a T5 channel problem. The extra hits in bars
3, 5 and 7 may be bar-6 hits recorded in its neighbours.

The data **x** positions reach ±10 cm, past the end of the longest bar
(±61.5 mm), and their spread (36.9 mm) is larger than any beam passing through
T5 can give (≲ 29 mm). Most likely these values include T5's x resolution
(roughly 2 cm), so fit the beam width to y or to bar fractions, not to x.

An earlier setting, σ = 3.2 cm in a 12 cm circle, gave an effective width of
only 2.55 cm. It is too narrow for this data.

## T5 position in the geometry

T5 sits inside the beam pipe, just upstream of the window: its downstream face
is 5 mm before the window. This placement comes from Martiloris/WCSim
(`a449968`). Before, it was at the tank-wall end of the pipe. A particle
started at the usual WCTE vertex (`/gps/position 0 -42.47625 -134.876 cm`,
direction +z) begins just downstream of the window, so T5 and window material
do not affect it in either placement. They only matter if the start point is
placed further upstream.

## Validation

All runs used `wcsim_runs/modified/WCTE.mac` geometry. Test files are in
`/lustre07/scratch/taghayor/wcsim_runs/beamspot_test/`.

- Single geantinos fired back through T5 at chosen offsets: bar length lies
  along global x, bars are stacked along +y, and each bar is hit as predicted.
- 5000 events, σ = 5 cm, T5 acceptance on: every start point lies inside a bar,
  z is unchanged, and the per-bar fractions are as in the table above.
- 2000 π⁺ at T = 198 MeV with a 2 % spread: mean p 307.41 MeV/c (expected
  307.37), σ_p 6.13 ± 0.10 MeV/c (expected 6.15), σ_T 5.58 MeV (expected 5.60),
  Gaussian shape.
- With every option unset, the start point and energy are identical to before.
- The WCSim truth output (`GetVtx`, beam energy) stores the spread values.
