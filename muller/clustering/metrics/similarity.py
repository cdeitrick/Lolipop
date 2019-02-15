import math
from typing import Union

import pandas
from dataclasses import dataclass


@dataclass
class PairCalculation:
	# Holds the values used to calculate the p-value for every pair.
	label: str
	pvalue: float
	X: float


def filter_out_invalid_timepoints(df: pandas.DataFrame, detected_cutoff: float, fixed_cutoff: float) -> pandas.DataFrame:
	"""
		Removes all timepoints from the dataframe where both values were either undetected or fixed.
	Parameters
	----------
	df:pandas.DataFrame
		A table with two series aligned along their timepoints.
	detected_cutoff: float
	fixed_cutoff: float

	Returns
	-------
	pandas.DataFrame
		A dataframe with any invalid timepoints removed.
	"""
	not_detected_fixed_df = df[df.lt(fixed_cutoff).any(axis = 1) & df.gt(detected_cutoff).any(axis = 1)]
	return not_detected_fixed_df


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
	not_detected_fixed_df = filter_out_invalid_timepoints(df, detected_cutoff, fixed_cutoff)
	# Remove timepoints where at least one trajectory was not fixed or was undetected.
	if not_detected_fixed_df.empty:
		p_value = calculate_overlap(left, right, fixed_cutoff)
		X = math.nan
	else:
		X = calculate_similarity(not_detected_fixed_df)
		p_value: float = 1.0 - math.erf(X)
	value = PairCalculation(
		label = point_label,
		X = X,
		pvalue = p_value
	)
	return value


def calculate_overlap(left: pandas.Series, right: pandas.Series, fixed_cutoff: float) -> float:
	"""
		Calculates the overlap of two series that are both undetected or fixed at all timepoints.
	Parameters
	----------
	left:pandas.Series
	right:pandas.Series
	fixed_cutoff: float

	Returns
	-------
	float
	"""
	left_fixed: pandas.Series = left[left.gt(fixed_cutoff)]
	right_fixed: pandas.Series = right[right.gt(fixed_cutoff)]

	if left_fixed.empty and right_fixed.empty:
		# Both are undetected, since they both failed the invalid timepoint filter.
		p_value = 1.0
	else:
		# Check if the trajectories overlap
		overlap = set(left_fixed.index) & set(right_fixed.index)
		# the p_value should be 1 if the fixed trajectories overlapped, and 0 otherwise.
		p_value = float((len(left_fixed) > 2 and len(right_fixed) > 2 and len(overlap) > 2))

	return p_value


def calculate_similarity(not_detected_fixed_df: pandas.DataFrame):
	# Find the mean frequency of each timepoint
	# index is timepoints,  values are frequencies
	not_detected_fixed_df = not_detected_fixed_df
	mean: pandas.Series = not_detected_fixed_df.mean(axis = 1)

	# Calculate sigma_freq
	# E(sigma) = (1/n) sum(sigma) = (1/n) sum(np(1-p)) == sum(p(1-p)
	# E(sigma_p) = (1/n) E(sigma) == 1/n(sum(p(1-p))
	# E(d_bar) = 1/n(sum(di)) == 1/n (n*sum(di))
	# pandas.Series.radd is slow for some reason. Use '-' operator instead.
	# noinspection PyTypeChecker
	sigma_freq: pandas.Series = mean.mul(1 - mean)
	# Difference of frequencies at each timepoint
	# difference: pandas.Series = not_detected_fixed_df.iloc[:, 0] - not_detected_fixed_df.iloc[:, 1]
	difference = not_detected_fixed_df.diff(axis = 1).iloc[:, 1]
	sigma_pair: float = sigma_freq.sum() / len(mean)
	# Sum of differences
	difference_mean: float = abs(difference).sum()

	X = difference_mean / (math.sqrt(2 * sigma_pair))

	return X


if __name__ == "__main__":
	pass
