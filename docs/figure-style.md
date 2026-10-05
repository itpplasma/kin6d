# kin6d plot and animation style

Scope: figures and movies produced by `app/*_demo.f90` through fortplot.
These rules affect presentation only. Each plotted band must still come from
the interval data written next to the image; the pixels never carry the bound.

## Canvas and encoding

- Use one canvas size per layout: 2x2 panels `figsize=[13.2, 8.8]`
  (1320x880 px), 2x1 or 1x2 panels `[12.8, 6.4]`. Give whole-pixel sizes; an
  even pixel count avoids encoder padding.
- 24 fps and at least 8 s (>= 192 frames) for a continuous evolution. The
  current 5 s clips at 12-16 fps are too short to read and look choppy. Hold
  the last state for about 1 s, either with repeated frames or a dedicated
  closing frame.
- Leave encoding to fortplot `save_animation`: H.264 yuv420p, BT.709 tagged,
  CRF 18, one keyframe per second, faststart (fortplot branch `mp4-quality`).
- Call `tight_layout()` on every subplot figure. With the fortplot fix it
  reserves room for the right-most x and top-most y tick labels.

## Text

- Titles state what a panel shows plus at most one changing quantity
  (`t = 3.00`). Move certification status and bookkeeping ("candidate",
  "1% budget met", "event derivative unavailable") to one `suptitle` or the
  gallery caption. Never use more than two `|`-separated fields in a title.
- Write symbols with mathtext or Unicode, not ASCII names: `$\omega_{pe} t$`,
  `$x/\lambda_D$`, `$v/v_{th}$`, `τ`, rather than `omega_pe t`, `lambdaD`,
  `v_th`. Mathtext and Greek letters were checked in fortplot raster output. The
  arrow glyph `→` renders as a box, so avoid it.
- Dimensionless axes: write the normalisation (`x / λ_D`) and drop `[1]`.
  Dimensional axes: give SI units in brackets.
- Print numbers in titles with a fixed format that cannot show
  denormals or Fortran exponents: `es9.2` becomes `1.2e-03` through a helper,
  and values below 1e-300 print as `0`. `2.0E-323` currently appears in
  `particles-vlasov`.

## Colour

- Sequential fields (`pcolormesh`, densities, f(x,v)): `cmap='viridis'` with
  fixed `vmin`/`vmax` for the whole movie, plus a `colorbar(label=...)`
  (fortplot at the pin rejects `colorbar` inside subplot grids, so until
  that is fixed state the range in the axis label or caption). Do
  not put the colour range in the title. For signed fields (perturbations,
  potentials) use a diverging map centred at zero, with symmetric limits.
- Categorical lines: use the Okabe-Ito palette via explicit RGB triples
  - blue `[0.000, 0.447, 0.698]`, vermillion `[0.835, 0.369, 0.000]`,
  - bluish green `[0.000, 0.620, 0.451]`, orange `[0.902, 0.624, 0.000]`,
  - reddish purple `[0.800, 0.475, 0.655]`, sky blue `[0.337, 0.706, 0.914]`,
  - black for exact/analytic references.

  Keep one colour per role across all demos: numerical = blue, reference or
  analytic = black dashed, second method/species = vermillion, bound = band in
  the colour of its curve.
- Do not encode a distinction by colour alone. Also vary the line style or
  marker (numerical solid, reference dashed, initial state dotted grey).
- When curves coincide (numerical over reference), draw the reference last as
  a thin dashed black line so both stay visible. Two coloured lines drawn on
  top of each other blend into a third colour.

## Interval bands

- `fill_between(x, lower, upper, color=<curve colour>, alpha=0.25)` drawn
  before the curve, with no dotted edge lines. This needs fortplot branch
  `mp4-quality`. Before that fix, fortplot ignored `alpha` in `fill_between`.
  Also set `xlim`/`ylim` explicitly on `pcolormesh` subplots: without
  them, the tight layout cannot see the mesh limits and clips the edge
  tick labels. Keep the band label in the
  legend ("numerical ± certified bound").
- When the band is narrower than about 2 px at the plotted scale (most demos
  today), show a second panel with `(numerical − reference)` and the band
  `±radius` on a scale that resolves it, as `plasma-oscillation` already does.
  State in the legend which bound is plotted.

## Axes and frames

- Fix `xlim`/`ylim` once for the whole movie from the global data/bound
  range (plus 5% padding), set before the first frame. Never autoscale per
  frame.
- Choose limits from the data. The occupation panel in `quantum-correlations`
  uses 125 orbital indices but only about 5 carry weight. Use a log y axis
  or zoom in.
- Geometric x–y panels (spheres, clouds, trajectories) need
  `axis('equal')`. Spheres currently render as ellipses.
- Put legends outside the data (`loc` in an empty corner) or below the axes.
  Several panels place a legend on top of particles or curves.
- Empirical CDFs use `step`, not scatter markers.
- Remove panels that stay constant for the whole movie. Replace them with one
  sentence in the caption.

## Line widths and markers

- Data lines 1.5–2 px, reference lines 1.2 px dashed, markers 4–6 px. For
  particle clouds with more than about 500 points, use small markers
  with alpha 0.5, which avoids moiré and GIF bloat (the `particles-*` GIFs are
  10 MB).
