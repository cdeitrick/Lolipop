import itertools
import math
from typing import Dict, Iterable, List, Tuple

import pandas


from muller.clustering.metrics import distance
from muller.widgets import get_valid_points



class CalculateMutationalDistances:
	def __init__(self, metric: str, detection_cutoff: float, fixed_cutoff: float):
		self.dlimit = detection_cutoff
		self.flimit = fixed_cutoff
		self.metric = metric

	def run(self, trajectories: pandas.DataFrame) -> Dict[Tuple[str, str], float]:
		pairwise_keys = self._generate_pairwise_combinations(trajectories.index)
		pair_array = dict()
		for left, right in pairwise_keys:
			left_trajectory = trajectories.loc[left]
			right_trajectory = trajectories.loc[right]

			# We only care about the timepoints such that `detection_cutoff` < f < `fixed_cutoff.
			# For now, lets require that both timepoints are detected and not yet fixed.
			left_reduced, right_reduced = get_valid_points(left_trajectory, right_trajectory, self.dlimit, self.flimit, inner = False)

			if left_reduced.empty or right_reduced.empty:
				distance_between_series = self.calculate_fixed_series_distance(left_trajectory, right_trajectory)
			else:
				distance_between_series = distance.calculate_distance(left_reduced, right_reduced, self.metric)

			pair_array[left, right] = pair_array[right, left] = distance_between_series

		return pair_array

	@staticmethod
	def _generate_pairwise_combinations(keys: Iterable[str]) -> List[Tuple[str, str]]:
		""" Generates a list of all possible pair of elements contained in `keys`"""

		# noinspection PyTypeChecker
		pair_combinations: Iterable[Tuple[str, str]] = itertools.combinations(keys, 2)
		return list(pair_combinations)

	@staticmethod
	def _get_maximum_value(items: Iterable[float]) -> float:
		""" Extracts the maximum values contained within an iterable."""
		return max(filter(lambda s: not math.isnan(s), items))

	def calculate_fixed_series_distance(self, left: pandas.Series, right: pandas.Series) -> float:
		"""
			Calculates the overlap of two series that are both undetected or fixed at all timepoints.
		Parameters
		----------
		left:pandas.Series
		right:pandas.Series
		"""
		left_fixed: pandas.Series = left[left.gt(self.flimit)].dropna()
		right_fixed: pandas.Series = right[right.gt(self.flimit)].dropna()
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
