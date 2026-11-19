# Clarifications to the frozen V5.3 documentation

Approved by the requester on 17 November 2026. These clarifications accompany, and do not replace, the byte-preserved original files in ESWA_V5_3_EVIDENCE.zip.

1. In ENDPOINT_MONOTONICITY_NOTE.md, the general timing feature is s=t_p. It equals t90 when p=0.9, the value fixed in the supplied implementation. The note's simultaneous general-p and t90 wording should be read with this distinction.
2. If B_tau(t) is defined as the untruncated integral of b_tau, it need not equal the initial speed v after stopping. The effective impulse for the stopped trajectory is min(v, B_tau(t)); this effective quantity equals v after stopping. The original theorem and implementation already use stopped trajectories. No theorem or computational result is changed.
3. The MATLAB header's phrase "5% feature-error" is inconsistent with its actual assignment eta=0.05. The setting is an absolute magnitude-error half-width of 0.05 m/s^2, as stated in the prompt, code assignment, README and article; it is not a percentage.
4. The standard mode uses 10 replications with 100 planning and 200 test draws. Its Student-t intervals use the source's rounded critical multiplier 2.262 for 9 degrees of freedom. Frozen interval endpoints are retained.

The two publication figures were redrawn solely from the unchanged source CSVs. Their deterministic plotting script is supplied separately. No MATLAB experiment, vehicle simulation or new random sampling was run during this revision.
