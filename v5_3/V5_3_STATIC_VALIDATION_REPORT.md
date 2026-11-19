# V5.3 static validation report

## Scope
This report records source-level checks only. MATLAB/Octave was not available
in the generation environment, so no claim of successful MATLAB execution
is made.

## Implemented checks and safeguards
- Function entry point matches the `.m` filename.
- Three V5.2 synthetic scenario parameter sets are retained.
- The fitted-proxy upper tau endpoint receives only the same sub-microsecond
  decimal-rounding correction as V5.2.
- The endpoint method checks `p >= 0.5`, nonnegative onset, nonempty supports,
  and ramp completion before stopping.
- The posterior interval update checks that each generated synthetic truth
  remains in the compatible interval.
- Outer-snapped tau intervals are checked to contain the continuous interval.
- The snapped-support robust objective is checked not to fall below the
  exact-support objective beyond numerical tolerance.
- The root-resolved time calculation evaluates phase boundaries and all
  real roots of piecewise quadratic relative velocity.
- Planning and test draws use separate seeds; the same seed protocol is
  applied to every tau-grid resolution.
- Monte Carlo intervals are labelled as within-scenario resampling
  variability, not real-fleet population intervals.
- All figures carry a synthetic-model footer.

## Outputs
- tau-grid support conservatism;
- decision-aware and fair random-expected policy convergence;
- action selection and STOP frequencies;
- time-grid convergence;
- factorised numerical error table.

## Limitations
- The validation remains within the fixed-cap/fixed-t90 ramp family.
- The three fleets are frozen synthetic designs.
- The fitted-proxy case combines configuration-specific fitted anchors
  from one vehicle and is not a simultaneous three-truck experiment.
- Root solving uses floating-point arithmetic, not interval arithmetic.
- No new measured data or field validation are introduced.
