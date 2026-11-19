"""Reduce the frozen continuous boxes under exact feature reports.

Each feasible shape has a known cap. An interior sample or t90 is one
affine equation in (onset, duration). The remaining optimization is a
corner, a timing segment, or a single point. Interval arithmetic encloses
coefficients that involve logarithms and exponentials.
"""

import math
from dataclasses import dataclass
from fractions import Fraction

from mpmath import iv


RAMP_TIMING = Fraction(9, 10)
FIXED_TIMES = tuple(map(Fraction, ("0.35", "0.65", "0.95")))


def decimal(value):
    return Fraction(str(value))


def interval(value):
    if isinstance(value, Fraction):
        return iv.mpf(value.numerator) / value.denominator
    return value


def timing_factor(shape):
    return RAMP_TIMING if shape == "ramp" else iv.ln(10)


@dataclass(frozen=True)
class Sample:
    time: Fraction
    elapsed: Fraction
    truth_shape: str

    @property
    def phase(self):
        if self.elapsed <= 0:
            return "zero"
        if self.truth_shape == "ramp" and self.elapsed >= 1:
            return "cap"
        return "interior"

    def normalized(self):
        if self.phase == "zero":
            return Fraction(0)
        if self.phase == "cap":
            return Fraction(1)
        if self.truth_shape == "ramp":
            return self.elapsed
        return 1 - iv.exp(-interval(self.elapsed))

    def equation(self, shape):
        if shape == self.truth_shape:
            slope = interval(self.elapsed)
        elif shape == "ramp":
            slope = self.normalized()
        else:
            slope = -iv.ln(1 - interval(self.normalized()))
        return interval(self.time), slope


@dataclass(frozen=True)
class ExactReport:
    cap: Fraction
    lag: Fraction
    duration: Fraction
    truth_shape: str
    timing: object
    samples: tuple[Sample, ...]

    def feature_values(self):
        features = [self.cap]
        if self.timing is not None:
            features.append(self.timing)
        for sample in self.samples:
            ratio = sample.normalized()
            features.append(self.cap * ratio if isinstance(ratio, Fraction) else interval(self.cap) * ratio)
        return features


def exact_report(truth, contract, extra_time):
    cap, lag, duration = map(decimal, truth[:3])
    shape = truth[3]
    timing = None
    if contract in {"D3", "D3+"}:
        timing = lag + RAMP_TIMING * duration if shape == "ramp" else interval(lag) + iv.ln(10) * interval(duration)
    times = FIXED_TIMES if contract == "D4" else ()
    if contract == "D3+" and extra_time is not None:
        times = (decimal(extra_time),)
    samples = tuple(Sample(time, (time - lag) / duration, shape) for time in times)
    return ExactReport(cap, lag, duration, shape, timing, samples)


def cell_index(value, quantum):
    """Floor an exact rational or certify that a transcendental avoids the boundary."""
    if isinstance(value, Fraction):
        return value // quantum
    index_interval = value / interval(quantum)
    index = math.floor(float(index_interval.mid))
    if not (index_interval.a >= index and index_interval.b < index + 1):
        raise ArithmeticError("Unresolved observation cell at available precision")
    return index


@dataclass(frozen=True)
class Extremum:
    shape: str
    lag: object
    duration: object
    method: str
    duration_range: tuple | None = None


def interval_max(left, right):
    return iv.mpf([max(left.a, right.a), max(left.b, right.b)])


def interval_min(left, right):
    return iv.mpf([min(left.a, right.a), min(left.b, right.b)])


def within(candidate, low, high):
    if candidate.b < low.a or candidate.a > high.b:
        return False
    if candidate.a >= low.b and candidate.b <= high.a:
        return True
    raise ArithmeticError("A candidate intersects a support boundary at available precision")


def exact_extrema(spec, report):
    """Return one maximizer for each feasible shape in the declared design."""
    lag_low, lag_high = map(decimal, spec["lag"])
    zero_times = [sample.time for sample in report.samples if sample.phase == "zero"]
    lag_low = max([lag_low, *zero_times])
    cap_times = [sample.time for sample in report.samples if sample.phase == "cap"]
    cap_time = min(cap_times) if cap_times else None
    interior = [sample for sample in report.samples if sample.phase == "interior"]
    assert decimal(spec["a"][0]) <= report.cap <= decimal(spec["a"][1])

    for shape in ("ramp", "first_order"):
        if lag_low > lag_high or (shape == "first_order" and cap_time is not None):
            continue
        if len(interior) == 3 and shape != report.truth_shape:
            # A straight ramp and a strictly concave exponential cannot agree
            # at three distinct times when all reported magnitudes are interior.
            continue

        duration_low, duration_high = map(decimal, spec[shape])
        lower_lag, upper_lag = map(interval, (lag_low, lag_high))
        lower_duration, upper_duration = map(interval, (duration_low, duration_high))
        equations = []
        if report.timing is not None:
            equations.append((interval(report.timing), interval(timing_factor(shape))))
        equations.extend(sample.equation(shape) for sample in interior)

        if not equations:
            # In this design any cap constraints in the no-equation case are
            # redundant on the prior box. A new design must check that premise.
            if cap_time is not None and lag_high + duration_high > cap_time:
                raise ValueError("Nonredundant plateau constraint on a two-dimensional support")
            yield Extremum(shape, upper_lag, upper_duration, "slow_corner", (lower_duration, upper_duration))
        elif len(equations) == 1:
            time, slope = equations[0]
            lower = interval_max(lower_duration, (time - upper_lag) / slope)
            upper = interval_min(upper_duration, (time - lower_lag) / slope)
            if cap_time is not None:
                assert shape == "ramp" and slope.b < 1
                upper = interval_min(upper, (interval(cap_time) - time) / (1 - slope))
            if lower.a > upper.b:
                continue
            if lower.b > upper.a:
                raise ArithmeticError("Unresolved timing-segment endpoint order")
            # Fixed timing orders cumulative braking along this segment. The
            # one-interior-sample ramp case also has slope >= 1/2 in this design.
            assert slope.a >= (0.5 if shape == "ramp" else 1)
            yield Extremum(shape, time - slope * lower, lower, "ordered_segment", (lower, upper))
        else:
            time_1, slope_1 = equations[0]
            time_2, slope_2 = equations[1]
            difference = slope_2 - slope_1
            assert difference.b < 0 or difference.a > 0
            if shape == report.truth_shape:
                yield Extremum(shape, interval(report.lag), interval(report.duration), "unique_truth")
                continue
            duration = (time_2 - time_1) / difference
            lag = time_1 - slope_1 * duration
            if not within(duration, lower_duration, upper_duration) or not within(lag, lower_lag, upper_lag):
                continue
            if cap_time is not None and not within(lag + duration, iv.mpf(0), interval(cap_time)):
                continue
            yield Extremum(shape, lag, duration, "unique_other_shape")
