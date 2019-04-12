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


def fixed_overlap(left: pandas.Series, right: pandas.Series, fixed_cutoff: float) -> float:
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
	left_fixed: pandas.Series = left[left.gt(fixed_cutoff)].dropna()
	right_fixed: pandas.Series = right[right.gt(fixed_cutoff)].dropna()

	if left_fixed.empty and right_fixed.empty:
		# Both are undetected, since they both failed the invalid timepoint filter.
		p_value = False
	else:
		# Since these genotypes fixed immediately, they should only be grouped together if they
		# fixed at the same timepoint.
		fixed_at_same_time = left_fixed.index[0] == right_fixed.index[0]
		# Check if the trajectories overlap completely
		overlap = set(left_fixed.index) & set(right_fixed.index)
		complete_overlap = len(overlap) == len(left_fixed) and len(overlap) == len(right_fixed)
		p_value = fixed_at_same_time and complete_overlap
	result = 0 if p_value else math.nan
	return result


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
		detected_points = pandas.DataFrame(
			{
				'left':  left_trajectory,
				'right': right_trajectory
			}
		)
		left_trajectory = detected_points['left']
		right_trajectory = detected_points['right']

		if left_trajectory.empty:
			distance_between_series = fixed_overlap(left_trajectory, right_trajectory, fixed_cutoff)
		else:
			distance_between_series = distance.calculate_distance(left_trajectory, right_trajectory, metric)
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
