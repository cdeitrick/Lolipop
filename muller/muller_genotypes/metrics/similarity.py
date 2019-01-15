import math
from typing import Optional, Union

import pandas
from dataclasses import dataclass


@dataclass
class PairCalculation:
	# Holds the values used to calculate the p-value for every pair.
	label: str
	pvalue: float
	X: float
	sigma: float
	difference_mean: float
	mean_series: Optional[pandas.Series]
	difference_series: Optional[pandas.Series]
	sigma_series: Optional[pandas.Series]


# noinspection PyTypeChecker
def calculate_p_value(left: pandas.Series, right: pandas.Series, detected_cutoff: float, fixed_cutoff: float) -> Union[float, PairCalculation]:
	"""
		Calculates the relative similarity between all trajectory pairs.

		these are consistent with n independent draws from a normal distribution,
        assuming an average variance of  n_{binom} p (1-p), where p is the
        average of the two frequencies, and n_{binom} is picked arbitrarily as a
        value that gives reasonable uncertainties given our data

        n random draws will have
        \sigma_{tot}^2 = \sum_i \sigma_i^2 =  n_{binom} \sum_i p_i(1 - p_i), from
        a property of normal distributions

        for \bar{X} = \sum_{i}^{n_X} X/n_X

        \sigma_{\bar{X}} = \sigma_{tot}/n_X = \sqrt( \sum_i p_i(1-p_i))/\sqrt{n_X}

        Finally, given \sigma and \bar{X} we can construct a p-value for
        the measurement by numerically computing the following integral:

        1 - \int_{- \bar{X}/ \sigma_{\bar{X}}}^{\bar{X}/
        \sigma_{\bar{X}}} \frac{1}{\sqrt{2 \pi}} e^{-x^2/2} dx
	Parameters
	----------
	left, right: pandas.Series
		- index: int
			Timepoints
		- values: float
			Frequencies
	detected_cutoff: float
	fixed_cutoff: float


	Returns
	-------
	"""
	# Merge into a dataframe for convienience
	if not isinstance(left, pandas.Series):
		left = pandas.Series(left)
	if not isinstance(right, pandas.Series):
		right = pandas.Series(right)
	point_label = f"{left.name} - {right.name}"
	df = pandas.concat([left, right], axis = 1)

	# Remove timepoints where at least one trajectory was not fixed or undetected.

	not_detected_fixed_df = df[df.lt(fixed_cutoff).any(axis = 1) & df.gt(detected_cutoff).any(axis = 1)]
	# Remove timepoints where at least one trajectory was not fixed or undetected.

	if not_detected_fixed_df.empty:
		left_fixed: pandas.Series = left[left.gt(fixed_cutoff)]
		right_fixed: pandas.Series = right[right.gt(fixed_cutoff)]

		if left_fixed.empty and right_fixed.empty:
			# Both are undetected
			p_value = 1.0
		else:
			# Check if the trajectories overlap
			overlap = set(left_fixed.index) & set(right_fixed.index)
			# the p_value should be 1 if the fixed trajectories overlapped, and 0 otherwise.
			p_value = float((len(left_fixed) > 2 and len(right_fixed) > 2 and len(overlap) > 2))

		value = PairCalculation(
			label = point_label,
			pvalue = p_value,
			X = 0 if p_value else math.inf,
			sigma = math.nan,
			difference_mean = math.nan,
			mean_series = None,
			difference_series = None,
			sigma_series = None
		)
	else:
		# Find the mean frequency of each timepoint
		# index is timepoints,  values are frequencies
		mean: pandas.Series = not_detected_fixed_df.mean(axis = 1)

		# Calculate sigma_freq
		# E(sigma) = (1/n) sum(sigma) = (1/n) sum(np(1-p)) == sum(p(1-p)
		# E(sigma_p) = (1/n) E(sigma) == 1/n(sum(p(1-p))
		# E(d_bar) = 1/n(sum(di)) == 1/n (n*sum(di))
		# pandas.Series.radd is slow for some reason. Use '-' operator instead.
		sigma_freq: pandas.Series = mean.mul(1 - mean)
		# Difference of frequencies at each timepoint
		# difference: pandas.Series = not_detected_fixed_df.iloc[:, 0] - not_detected_fixed_df.iloc[:, 1]
		difference = not_detected_fixed_df.diff(axis = 1).iloc[:, 1]
		sigma_pair: float = sigma_freq.sum() / len(mean)
		# Sum of differences
		difference_mean: float = abs(difference).sum()

		X = difference_mean / (math.sqrt(2 * sigma_pair))

		p_value: float = 1.0 - math.erf(X)

		value = PairCalculation(
			label = point_label,
			pvalue = p_value,
			X = X,
			sigma = sigma_pair,
			difference_mean = difference_mean,
			mean_series = mean,
			difference_series = difference,
			sigma_series = sigma_freq
		)

	return value


if __name__ == "__main__":
	pass
