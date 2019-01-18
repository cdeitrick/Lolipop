import pandas
from typing import Dict, Tuple, List
import itertools
try:
	from muller_genotypes.metrics.similarity import PairCalculation, calculate_p_value
	from muller_genotypes.metrics.dynamic_time_warping import calculate_dtw
except ModuleNotFoundError:
	from .similarity import PairCalculation, calculate_p_value
	from .dynamic_time_warping import calculate_dtw

def calculate_pairwise_metric(trajectories: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float, metric: str) -> Dict[Tuple[str,str], PairCalculation]:
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
		if metric == "similarity":
			calculation = calculate_p_value(left_trajectory, right_trajectory, detection_cutoff, fixed_cutoff)
		elif metric == 'dtw':
			calculation, *_ = calculate_dtw(left_trajectory, right_trajectory)
			calculation = PairCalculation(label = f"{left_trajectory.name} - {right_trajectory.name}", pvalue = None, X = calculation)
		else:
			message = f"'{metric}' is not an available metric."
			raise ValueError(message)
		pair_array[left, right] = pair_array[right, left] = calculation

	return pair_array