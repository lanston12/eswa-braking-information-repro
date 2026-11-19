# V5.3 analytical endpoint note

For the ramp response
b_tau(t) = a * min(max((t - s + p*tau)/tau, 0), 1),
with fixed a, s=t90 and p in [0.5,1], define cumulative braking
B_tau(t) = integral_0^t b_tau(u) du.

During the ramp, set q=t-s. Then
B_tau(t) = a*(q+p*tau)^2/(2*tau), and
dB_tau/dtau = a*(q+p*tau)*(p*tau-q)/(2*tau^2).

On the ramp interval, q+p*tau >= 0 and q <= (1-p)*tau, hence
p*tau-q >= (2p-1)*tau >= 0.

During the plateau,
B_tau(t) = a*(t-s+(p-1/2)*tau), so
dB_tau/dtau = a*(p-1/2) >= 0.

Before onset B_tau=0; after stopping B_tau=v. Therefore cumulative
braking is nondecreasing in tau. The corresponding velocity and
displacement are nonincreasing in tau.

For independent product supports, the worst directed closure
max_t[x_follower(t)-x_leader(t)] is therefore attained by the
largest admissible tau for the leader and the smallest admissible
tau for the follower. This justifies endpoint evaluation within
the stated ramp family and assumptions.

This is a model-specific analytical monotonicity argument. It does
not apply automatically to first-order, arbitrary measured, or
correlated response supports.
