# Decision-relevant braking information: minimal reproduction materials

This repository accompanies an ESWA submission on braking information and robust truck-platoon formation. It contains a small, inspectable subset of the author-held reproduction archives. It does not contain a new vehicle experiment, an operational safety certificate, or the complete original archives.

## Quick verification (Python standard library)

From the repository root, run:

```sh
python verify_reproduction.py
```

The script checks source SHA-256 values, the exact 6.1 m equal-t90 witness, the Figure 3 measurement inversion, finite-study row counts, V4 replication-to-summary means and the primary 0.1663 m-equivalent paired advantage, the separate V5.3 summaries, empirical event-to-summary consistency, mixed empirical outcomes, and the missing Cambridge 18b endpoint. It writes `verification_report.json`. These checks reconcile frozen results; they do not rerun all original simulations.

## Recreate the four Supplement N figures (MATLAB)

In MATLAB, change directory to `empirical` and run:

```matlab
ESWA_Empirical_Evidence_Enhancement_OneClick
```

The script reads the four supplied, full-precision frozen CSV files. It produces four PDF/PNG/FIG figures, four plotting CSV files and a scope note in `ESWA_EMPIRICAL_EVIDENCE_OUTPUT`. No model is fitted again. The four plots cover event-level paired errors, pressure-reference threshold sensitivity, near-optimal t90 profile widths, and retrospective Cambridge prefix extrapolation. The script was executed using MATLAB R2022b on Windows; see `VALIDATION.md` for the actual verification scope. Generated files are ignored by Git.

## Original simulation source and frozen results

- `synthetic`: the frozen finite-support design, physics and planner code, exact certificate and continuous-enclosure code, and selected frozen results. The 576 policy rows comprise 24 synthetic fleets; 2,100 cost rows belong to a separate cost sensitivity calculation. Dependencies in the original `requirements.txt` and environment record are preserved as historical records, not a freshly verified installation recipe. For a new finite-study run, install NumPy and the dependencies imported by the selected script, then run `python synthetic/code/run_targeted.py --workers 2 --output synthetic/recomputed` from the repository root. A complete run can be costly. This command's full-study execution is NOT VERIFIED in this release.
- `v4`: the unchanged one-query MATLAB source and frozen CSVs. The primary overlap result uses 30 paired replications. Objective gains are **m-equivalent**, not measured physical distances. Both control designs select STOP, and forced variance-first queries can lose value.
- `v5_3`: the unchanged, separate ten-replication resolution study and its source CSVs. Read `DOCUMENTATION_CLARIFICATIONS.md` alongside the original frozen note; the clarification distinguishes raw and speed-truncated braking impulse and general t_p from t90. These results do not replace the primary V4 experiment.
- `explanatory_figures`: editable TikZ sources. These illustrate explicitly stated models and are not measured observations.

The original V4 and V5.3 MATLAB programs may run substantial new calculations. Their complete Monte Carlo reruns are NOT VERIFIED in this release. Run each in its own scratch directory to preserve the supplied frozen evidence. The full author archives, including the complete historical output and raw-record decoding/fitting workflow, are available on reasonable request from the corresponding author, Jianbo Feng, at fengjianbo@bucea.edu.cn.

## Data provenance and limits

Public vehicle recordings originate from the [Colorado State University correlated DDEC VI dataset](https://www.engr.colostate.edu/~jdaily/tucrrc/CorrelatedDDEC6DataSet.html) and [Henderson and Cebon's Cambridge dataset](https://doi.org/10.17863/CAM.69067). The repository redistributes derived frozen CSV results, not raw third-party recordings. Use the original providers for source files and their terms. Cambridge publication-selected traces and within-record prefix extrapolation are not independent fleet validation. The reported empirical mixed results, failed endpoints and synthetic negative controls are intentionally retained.

`PROVENANCE.json` records the source and SHA-256 of every carried source/evidence file. `SHA256SUMS.txt` also covers the release's support scripts and documentation. Frozen originals have not been numerically edited. A separate authoring revision reads exact CSV values instead of rounded literals in the new empirical plotting script.

OpenAI Codex (GPT-6.1 Sol, as confirmed by the authors) assisted code and explanatory-figure development and the assembly of this package. The authors reviewed and edited AI-assisted output as needed and accept responsibility for the article. No additional software licence is granted by this release; third-party materials retain their original terms.

This public repository identifies its owner. The submission also has a separate anonymous manuscript in which this repository link is withheld for review.
