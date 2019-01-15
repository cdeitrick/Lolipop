import itertools
from typing import Dict, List, Tuple, Union

import pandas

try:
	from muller_genotypes.metrics.similarity import calculate_p_value, PairCalculation
except ModuleNotFoundError:
	from .similarity import calculate_p_value, PairCalculation

PairwiseArrayType = Dict[Tuple[str, str], PairCalculation]


class PairwiseCalculation:
	"""
		Calculates and holds the calculations for all pairwise elements.

		Usage
		-----
		calculation = PairwiseCalculation()
		# Update values
		pair_array = calculation.update_values(timeseries, detection_cutoff, fixed_cutoff)

	"""
	def __init__(self):
		self.pairwise_values: PairwiseArrayType = None

	def asitem(self, itemtype:str)->Dict[Tuple[str,str], float]:
		""" Returns a pairwise dictionary mapping all key to the attribute in PairCalculation defined by `itemtype`."""
		return {(left, right):self.get(left, right ,itemtype = itemtype) for left, right in self.pairwise_values}

	def update_values(self, timeseries, detection_cutoff: float, fixed_cutoff: float, metric: str = 'similarity'):
		"""
			Generates a dictionary mapping every 2-pair combination of trajectories to their corresponding distance/probability calculation.
			This dictionary is cached by the parent `PairwiseCalculation` object.
		Parameters
		----------
		timeseries: pandas.DataFrame
			A table of trajectories, where the index is the trajectory label and columns indicate timepoints.
		detection_cutoff: float
		fixed_cutoff: float
		metric: str; default 'similarity'
			The distance metric to use. Currently only one metric is implemented.

		Returns
		-------
		PairwiseArrayType
		"""
		if self.pairwise_values:
			_current_trajectory_labels = set(timeseries.index)
			pair_array = {
				k: v for k, v in self.pairwise_values.items()
				if (k[0] in _current_trajectory_labels and k[1] in _current_trajectory_labels)
			}
			self.pairwise_values = pair_array.copy()
		else:
			pair_array = calculate_pairwise_metric(
				timeseries,
				detection_cutoff = detection_cutoff,
				fixed_cutoff = fixed_cutoff,
				metric = metric
			)
			self.pairwise_values = pair_array.copy()
		return pair_array

	def squareform(self, itemtype:str):
		""" Converts a dictionary with all pairwise values for a set of points into a square matrix representation.

		Parameters
		----------
		itemtype: str
			One of the available PairCalculation attributes.
		"""
		keys = sorted(set(itertools.chain.from_iterable(self.pairwise_values.keys())))
		_square_map = dict()
		for left in keys:
			series = dict()
			for right in keys:
				value = self.get(left, right, itemtype, 0)
				series[right] = value
			_square_map[left] = series
		return pandas.DataFrame(_square_map)

	def get(self, left, right, itemtype:str = None, default = None)->Union[PairCalculation, float]:
		result = self.pairwise_values.get((left, right), default)
		if itemtype:
			result = getattr(result, itemtype)
		return result

def calculate_pairwise_metric(trajectories: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float, metric: str) -> PairwiseArrayType:
	"""
	Parameters
	----------
	trajectories: pandas.DataFrame
		A table of mutational trajectories. Should be a normal trajectory table.
	detection_cutoff: float
	fixed_cutoff: float
	metric: str

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
		if metric == "similarity":
			calculation = calculate_p_value(left_trajectory, right_trajectory, detection_cutoff, fixed_cutoff)
		else:
			message = f"'{metric}' is not an available metric."
			raise ValueError(message)
		pair_array[left, right] = pair_array[right, left] = calculation

	return pair_array
