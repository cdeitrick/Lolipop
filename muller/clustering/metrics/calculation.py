import itertools
import math
from typing import Dict, List, Tuple

import pandas
from loguru import logger

try:
	from muller.clustering.metrics import distance
	from muller import widgets
except ModuleNotFoundError:
	from . import distance
	from ... import widgets
try:
	from tqdm import tqdm
except ModuleNotFoundError:
	tqdm = None

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
	if left_fixed.empty or right_fixed.empty:
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
	metric: {'similarity', 'binomial', 'pearson', 'minkowski'}

	Returns
	-------
	dict of PairArrayValue
	Each key in the dictionary corresponds to a pair of trajectory ids which map to the p-value for that pair.
	The order of ids does not matter.
	"""
	logger.debug("Calculating the pairwise values...")
	logger.debug(f"\t detection limit: {detection_cutoff}")
	logger.debug(f"\t fixed limit: {fixed_cutoff}")
	logger.debug(f"\t metric: {metric}")

	# noinspection PyTypeChecker
	pair_combinations: List[Tuple[str, str]] = list(itertools.combinations(trajectories.index, 2))
	pair_array = dict()
	if len(pair_combinations) > 10000 and tqdm:
		progress_bar = tqdm(total = len(pair_combinations))
	else:
		progress_bar = None
	for left, right in pair_combinations:
		if progress_bar:
			progress_bar.update(1)
		left_trajectory = trajectories.loc[left]
		right_trajectory = trajectories.loc[right]

		# We only care about the timepoints such that `detection_cutoff` < f < `fixed_cutoff.
		# For now, lets require that both timepoints are detected and not yet fixed.
		# There is an issue related to comparing fixed genotypes against non-fixed genotypes.
		# These were grouped together:
		# recG    0    0.14    0.319    1    1    1    1
		# PA14_RS20565< 0    0.153    0.231    0    0    0    0
		#left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, detection_cutoff, fixed_cutoff, inner = False)
		left_was_fixed = widgets.fixed(left_trajectory, fixed_cutoff)
		right_was_fixed = widgets.fixed(right_trajectory, fixed_cutoff)
		if left_was_fixed == right_was_fixed:
			left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, detection_cutoff, fixed_cutoff, inner = False)
		else:
			left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, detection_cutoff, inner = False)

		if left_reduced.empty or right_reduced.empty:
			# Treat both trajectories as fixed immediately.
			distance_between_series = fixed_overlap(left_trajectory, right_trajectory, fixed_cutoff)
		else:
			distance_between_series = distance.calculate_distance(left_reduced, right_reduced, metric)

		pair_array[left, right] = pair_array[right, left] = distance_between_series
	if progress_bar: progress_bar.close()
	# Assume that any pair with NAN values are the maximum possible distance from each other.

	maximum_distance = max(filter(lambda s: not math.isnan(s), pair_array.values()))
	pair_array = {k: (v if not math.isnan(v) else maximum_distance) for k, v in pair_array.items()}

	return pair_array
