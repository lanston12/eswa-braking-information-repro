# Release validation (19 November 2026)

- PASS: all 13 checks in `verify_reproduction.py`, including 60 carried-file SHA-256 values, exact witness arithmetic, frozen means, negative controls and empirical mixed outcomes.
- PASS: MATLAB R2022b actually executed the empirical script from this package directory with process exit code 0; all four FIG files were reopened successfully. Four PDF/PNG/FIG and four CSV outputs were generated.
- PASS: synthetic planner imports and the 24-fleet frozen design loads with its original design hashes.
- NOT VERIFIED in this release: a new complete finite-study run; new V4 or V5.3 Monte Carlo simulations; reacquisition and decoding of all public raw recordings; MATLAB versions other than R2022b.

The paper's primary V4 and separate V5.3 results are reconciled against their frozen individual replications. No frozen experimental values were improved or replaced. Verification logs and complete source archives are retained by the authors.
