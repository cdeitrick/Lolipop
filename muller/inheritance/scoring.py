
import pandas
import math
from loguru import logger
from widgets import get_valid_points

try:
	from inheritance import order_by_area
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


	score = int(difference.sum() > (cutoff * len(difference)/2))
	logger.debug(f"{nested_genotype.name}\t{unnested_genotype.name}\t{difference.sum():.2f}\t{cutoff * len(difference)/2:.2f}")
	# There may be a point where the unnested genotype is greater than the nested genotype.

	if any(difference < 0):
		score -= 0.5

	return score

def calculate_derivative_score(left: pandas.Series, right: pandas.Series, detection_cutoff:float, cutoff:float) -> float:
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
	df = get_valid_points(left, right, detection_cutoff, flimit = 0.97)
	if df.empty:
		covariance = math.nan
	else:
		covariance = df.iloc[:, 0].cov(df.iloc[:, 1])
	if math.isnan(covariance): score = 0
	elif covariance > cutoff:score = 2
	elif covariance < -cutoff:score = -2
	else: score = 0

	return score

def calculate_area_score(nested_genotype:pandas.Series, unnested_genotype:pandas.Series)->float:
	"""
		Calculates a score based on the probability the unnested genotype is a subset of the nested_genotype.
		This take into account both the area of the nested genotype and the area of all other genotypes not in the nested genotype.
	Parameters
	----------
	nested_genotype, unnested_genotype: pandas.Series

	"""

	# noinspection PyTypeChecker
	other_genotypes:pandas.Series = 1 - nested_genotype

	common_area_nested = order_by_area.calculate_common_area(nested_genotype, unnested_genotype)
	common_area_other = order_by_area.calculate_common_area(other_genotypes, unnested_genotype)

	is_subset_nested = order_by_area.is_subset(nested_genotype, unnested_genotype)
	is_subset_other = order_by_area.is_subset(other_genotypes, unnested_genotype)

	if is_subset_nested and is_subset_other:
		score = int(common_area_nested > 2*common_area_other)
	elif is_subset_nested:
		score = 2
	else:
		score = -2

	return score




