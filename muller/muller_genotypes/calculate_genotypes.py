from typing import Dict, Optional, Tuple

import pandas
from dataclasses import dataclass

try:
	from muller_genotypes.metrics.similarity import PairwiseArrayType, PairCalculation, calculate_pairwise_trajectory_similarity
	from muller.muller_genotypes.methods import matlab_method, hierarchical_method
	from muller.muller_genotypes.genotype_average import calculate_mean_genotype
except ModuleNotFoundError:
	from .metrics.similarity import PairwiseArrayType, PairCalculation, calculate_pairwise_trajectory_similarity
	from .methods import matlab_method, hierarchical_method
	from .genotype_average import calculate_mean_genotype

PAIRWISE_P_VALUES: PairwiseArrayType = None
PAIRWISE_CALCULATIONS: Dict[Tuple[str, str], PairCalculation] = None
REMOVED_P_VALUES: PairwiseArrayType = dict()


@dataclass
class GenotypeOptions:
	detection_breakpoint: float  # Minimum frequency to be considered detected.
	fixed_breakpoint: float  # Frequency at which a mutation is considered fixed.
	n_binom: Optional[int]
	similarity_breakpoint: float  # The cutoff indicating two trajectories are related.
	# The cutoff indicating two trajectories that were originally sorted into the same genotype are not
	# actually related.
	difference_breakpoint: float
	method: str

	@classmethod
	def from_matlab(cls) -> 'GenotypeOptions':
		return GenotypeOptions(
			detection_breakpoint = 0.03,
			fixed_breakpoint = 0.97,
			n_binom = 5,
			similarity_breakpoint = 0.05,
			difference_breakpoint = 0.10,
			method = 'matlab'
		)


def generate_pair_array(timeseries: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float) -> PairwiseArrayType:
	# This is a global variable so that the p-values do not need to be re-computed for every iteration of the genotype filters.
	global PAIRWISE_P_VALUES
	global PAIRWISE_CALCULATIONS
	if PAIRWISE_P_VALUES:
		_current_trajectory_labels = set(timeseries.index)
		pair_array = {k: v for k, v in PAIRWISE_P_VALUES.items() if (k[0] in _current_trajectory_labels and k[1] in _current_trajectory_labels)}

	else:
		pair_array = calculate_pairwise_trajectory_similarity(
			timeseries,
			detection_cutoff = detection_cutoff,
			fixed_cutoff = fixed_cutoff
		)
		if PAIRWISE_CALCULATIONS is None:
			PAIRWISE_CALCULATIONS = {k: v for k, v in pair_array.items() if not isinstance(v, float)}

		pair_array = {k: v for k, v in pair_array.items()}

		PAIRWISE_P_VALUES = pair_array

	return pair_array


def workflow(timepoints: pandas.DataFrame, options: GenotypeOptions) -> Tuple[pandas.DataFrame, pandas.Series]:
	"""

	Parameters
	----------
	timepoints: pandas.DataFrame
		A timeseries dataframe, usually generated from `import_table.import_trajectory_table`.
			- Index -> str
				Names unique to each trajectory.
			- Columns -> int
				The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	options: GenotypeOptions
		An instance of `GenotypeOptions` with the desired options. Can be generated automatically via `GenotypeOptions.from_breakpoints(0.03)`.

	Returns
	-------
	pandas.DataFrame, pandas.Series
		- The genotype table
		- A map of genotypes to members.
	"""
	# calculate the similarity between all pairs of trajectories in the population.
	pair_array = generate_pair_array(timepoints, options.detection_breakpoint, options.fixed_breakpoint)

	if options.method == "matlab":
		genotypes = matlab_method(timepoints, pair_array, options.similarity_breakpoint, options.difference_breakpoint)
	elif options.method == "hierarchy":
		genotypes = hierarchical_method(pair_array, options.similarity_breakpoint)
	else:
		raise ValueError(f"Invalid clustering method: {options.method}")

	_mean_genotypes = calculate_mean_genotype(genotypes, timepoints)
	genotype_members = _mean_genotypes.pop('members')
	_mean_genotypes = _mean_genotypes[sorted(_mean_genotypes.columns)]

	return _mean_genotypes, genotype_members


if __name__ == "__main__":
	pass
