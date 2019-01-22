import itertools
from typing import Dict, List, Tuple

import pandas

try:
	from muller_genotypes.metrics import distance
except ModuleNotFoundError:
	from . import distance


def filter_out_invalid_timepoints(left: pandas.Series, right: pandas.Series, detected_cutoff: float, fixed_cutoff: float) -> Tuple[
	pandas.Series, pandas.Series]:
	"""
		Removes all timepoints from the dataframe where both values were either undetected or fixed.
	Parameters
	----------
	left:pandas.DataFrame
	right:pandas.DataFrame
	detected_cutoff:float
	fixed_cutoff:float

	Returns
	-------
	pandas.DataFrame
		A dataframe with any invalid timepoints removed.
	"""
	df = pandas.concat([left, right], axis = 1)
	not_detected_fixed_df = df[df.lt(fixed_cutoff).any(axis = 1) & df.gt(detected_cutoff).any(axis = 1)]
	return not_detected_fixed_df.iloc[:, 0], not_detected_fixed_df.iloc[:, 1]


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
	pair_combinations: List[Tuple[str, str]] = itertools.combinations(trajectories.index, 2)
	pair_array = dict()
	for left, right in pair_combinations:
		left_trajectory = trajectories.loc[left]
		right_trajectory = trajectories.loc[right]
		# point_label = f"{left_trajectory.name} - {right_trajectory.name}"
		left_trajectory, right_trajectory = filter_out_invalid_timepoints(left_trajectory, right_trajectory, detection_cutoff, fixed_cutoff)
		if left_trajectory.empty:
			calculation = calculate_overlap(left_trajectory, right_trajectory, fixed_cutoff)
		elif metric == "similarity":
			calculation = distance.binomial_probability(left_trajectory, right_trajectory)
		elif metric == 'pearson':
			calculation = distance.pearson_correlation_distance(left_trajectory, right_trajectory)
		elif metric == 'minkowski':
			calculation = distance.minkowski_distance(left_trajectory, right_trajectory, 2)
		elif metric == 'binomial':
			calculation = distance.binomial_distance(left_trajectory, right_trajectory)
		elif metric == 'binomialp':
			calculation = distance.binomial_probability(left_trajectory, right_trajectory)
		elif metric == 'dtw':
			calculation= distance.dynamic_time_warping(left_trajectory, right_trajectory)
		else:
			message = f"'{metric}' is not an available metric."
			raise ValueError(message)
		pair_array[left, right] = pair_array[right, left] = calculation

	return pair_array
