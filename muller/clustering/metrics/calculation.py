import itertools
import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__file__)
import pandas
import math

try:
	from clustering.metrics import distance
	from widgets import get_valid_points
except ModuleNotFoundError:
	from . import distance
	from ...widgets import get_valid_points


def normalize(series: pandas.Series) -> pandas.Series:
	mean = series.mean()
	sigma = series.var()
	if math.isnan(sigma) or len(series) < 2:
		normalized_series = None
	else:
		normalized_series = (series - mean) / sigma
	return normalized_series


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


def calculate_pairwise_metric(trajectories: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float, metric: str) -> Dict[
	Tuple[str, str], float]:
	"""
	Parameters
	----------
	trajectories: pandas.DataFrame
		A table of mutational trajectories. Should be a normal trajectory table.
	detection_cutoff: float
	fixed_cutoff: float
	metric: {'similarity', 'dtw'}

	Returns
	-------
	dict of PairArrayValue
	Each key in the dictionary corresponds to a pair of trajectory ids which map to the p-value for that pair.
	The order of ids does not matter.
	"""
	logger.info("Calculating the pairwise values...")
	logger.info(f"\t detection limit: {detection_cutoff}")
	logger.info(f"\t fixed limit: {fixed_cutoff}")
	logger.info(f"\t metric: {metric}")

	pair_combinations: List[Tuple[str, str]] = itertools.combinations(trajectories.index, 2)
	pair_array = dict()
	for left, right in pair_combinations:
		left_trajectory = trajectories.loc[left]
		right_trajectory = trajectories.loc[right]
		detected_points = get_valid_points(left_trajectory, right_trajectory, detection_cutoff, inner = False)
		#detected_points = pandas.concat([left_trajectory, right_trajectory], axis = 1)
		detected_points.columns = ['left', 'right']
		left_trajectory = detected_points['left']
		right_trajectory = detected_points['right']

		if left_trajectory.empty:
			distance_between_series = calculate_overlap(left_trajectory, right_trajectory, fixed_cutoff)
		elif metric == "similarity":
			distance_between_series = distance.binomial_probability(left_trajectory, right_trajectory)
		elif metric == 'pearson':
			distance_between_series = distance.pearson_correlation_distance(left_trajectory, right_trajectory)
		elif metric == 'minkowski':
			distance_between_series = distance.minkowski_distance(left_trajectory, right_trajectory, 2)
		elif metric == 'binomial':
			distance_between_series = distance.binomial_distance(left_trajectory, right_trajectory)
		elif metric == 'binomialp':
			distance_between_series = distance.binomial_probability(left_trajectory, right_trajectory)
		elif metric == 'jaccard':
			distance_between_series = distance.jaccard_distance(left_trajectory, right_trajectory)
		elif metric == "combined":
			distance_between_series_pearson = distance.pearson_correlation_distance(left_trajectory, right_trajectory)
			distance_between_series_minkowski = distance.minkowski_distance(left_trajectory, right_trajectory, 2)
			#distance_between_series = (distance_between_series_jaccard / 2) + distance_between_series_pearson
			distance_between_series = (2*distance_between_series_pearson) + distance_between_series_minkowski
		else:
			message = f"'{metric}' is not an available metric."
			raise ValueError(message)
		pair_array[left, right] = pair_array[right, left] = distance_between_series

	# Assume that any pair with NAN values are the maximum possible distance from each other.
	maximum_distance = max(filter(lambda s: not math.isnan(s), pair_array.values()))
	pair_array = {k: (v if not math.isnan(v) else maximum_distance) for k, v in pair_array.items()}
	# Log the values for debugging
	keys = sorted(set(i[0] for i in pair_array.keys()))
	for key in keys:
		items = {k: v for k, v in pair_array.items() if k[0] == key}
		for (l, r), v in sorted(items.items(), key = lambda s: s[-1], reverse = False):
			line = f"{l}\t{r}\t{v:.2E}"
			logger.debug(line)
	return pair_array
