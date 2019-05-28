import math

import pandas

from muller.widgets import get_valid_points

try:
	from muller.inheritance import order_by_area
except ModuleNotFoundError:
	from . import order_by_area


def calculate_additive_score(nested_genotype: pandas.Series, unnested_genotype: pandas.Series, cutoff: float) -> float:
	"""
		Tests whether the nested genotype is consistenty larger than  the nested genotype. The scoring is as follows:
		`nested_genotype` always larger than `unnested_genotype`: 1
		`nested_genotype` consistently larger than `unnested_genotype`: 0.5
		`nested_genotype` not consistently larger than `unnested_genotype`: 0
		`nested_genotype` sometimes smaller than `unnested_genotype`: -0.5
	Parameters
	----------
	nested_genotype, unnested_genotype: pandas.Series
	cutoff: The minimum distance between each timepoint between `nested_genotype` and `unnested_genotype`

	"""
	difference = nested_genotype - unnested_genotype

	score = int(difference.sum() > (cutoff * len(difference) / 2))

	# There may be a point where the unnested genotype is greater than the nested genotype.

	if any(difference < 0):
		score -= 1
	# logger.debug(f"{unnested_genotype.name}\t{nested_genotype.name}\t{score}")
	return score


def calculate_derivative_score(left: pandas.Series, right: pandas.Series, detection_cutoff: float, cutoff: float) -> float:
	"""
		Tests whther the two series are correlated or anticorrelated with each other. The scoring is as follows:
		correlated: 2
		uncorrelated: 0
		anticorrelated: -2
	Parameters
	----------
	left, right: pandas.Series
		The two series to test.
	detection_cutoff: float
		Used to determine which points to use for the comparison.
	cutoff: float
		The cutoff value to determine whther the series are actually correlated/anticorrelated. Any value within the range [-`cutoff`, `cutoff`]
		results in a score of 0.
	"""
	# Pandas implementation of the derivative check, since it basically just checks for covariance.
	valid_left, valid_right = get_valid_points(left, right, detection_cutoff, flimit = 0.97)
	if valid_left.empty:
		covariance = math.nan
	else:
		covariance = valid_left.cov(valid_right)

	if covariance > cutoff: score = 2
	elif covariance < -cutoff: score = -2
	else: score = 0
	from loguru import logger
	logger.debug(f"{right.name}\t{left.name}\t{score}\t{covariance}")
	return score


def calculate_area_score(nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
	"""
		Calculates a score based on the probability the unnested genotype is a subset of the nested_genotype.
		This take into account both the area of the nested genotype and the area of all other genotypes not in the nested genotype.
	Parameters
	----------
	nested_genotype, unnested_genotype: pandas.Series

	"""

	# If the nested genotype is not fixed, group the remaining frequencies into an `other` category.
	# noinspection PyTypeChecker
	other_genotypes: pandas.Series = 1 - nested_genotype
	area_nested = order_by_area.area_of_series(nested_genotype)
	area_unnested = order_by_area.area_of_series(unnested_genotype)
	area_other = order_by_area.area_of_series(other_genotypes)
	common_area_nested = order_by_area.calculate_common_area(nested_genotype, unnested_genotype)
	common_area_other = order_by_area.calculate_common_area(other_genotypes, unnested_genotype)

	is_subset_nested = order_by_area.is_subset(area_nested, area_unnested, common_area_nested)
	is_subset_other = order_by_area.is_subset(area_other, area_unnested, common_area_other)

	if is_subset_nested and is_subset_other:
		score = int(common_area_nested > 2 * common_area_other)
	elif is_subset_nested:
		score = 2
	else:
		score = -2
	# logger.debug(f"{unnested_genotype.name}\t{nested_genotype.name}\t{score}")
	return score


def calculate_subtractive_score(left: pandas.Series, right: pandas.Series, fixed_cutoff: float, cutoff: float) -> int:
	"""
		Tests whether two genotypes consistently sum to a value greater than the fixed breakpoint. This suggests that one of the genotypes
		is in the background of the other, since otherwise the maximum combined frequency should, at most, be equal to the fixed cutoff value.
	Parameters
	----------
	left, right: pandas.Series
	fixed_cutoff: The fixed breakpoint.
	cutoff: Governs whether two series consistently sum to a value greater than `fixed_cutoff`. Should be in the range [0,1]
	"""

	combined = left + right
	above_fixed = combined - fixed_cutoff

	# Finally, check whether the frequency of this series is consistently above the `cutoff` value.
	# Timepoints that do not sum to greater than the fixed breakpoint will be negative.
	result = above_fixed.mean() > cutoff

	return int(result)
