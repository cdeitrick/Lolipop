from typing import Any, Tuple, List

import pandas
from options import GenotypeOptions
try:
	from muller.muller_genotypes.average import calculate_mean_genotype
	from muller_genotypes.metrics import PairwiseCalculation, calculate_pairwise_metric
	from muller_genotypes.methods import calculate_genotypes_from_given_method
	from muller_genotypes.filters import filter_genotypes
	from muller.widgets import map_trajectories_to_genotype
except ModuleNotFoundError:
	from .average import calculate_mean_genotype
	from .metrics import PairwiseCalculation, calculate_pairwise_metric
	from .methods import calculate_genotypes_from_given_method
	from .filters import filter_genotypes
	from widgets import map_trajectories_to_genotype

PAIRWISE_CALCULATIONS = PairwiseCalculation()


def _update_pairwise_array(timepoints: pandas.DataFrame, options: GenotypeOptions):
	global PAIRWISE_CALCULATIONS
	if PAIRWISE_CALCULATIONS:
		# PAIRWISE_CALCULATIONS already contains previous calcuations.
		PAIRWISE_CALCULATIONS = PAIRWISE_CALCULATIONS.reduce(timepoints.index)
	else:
		pair_array = calculate_pairwise_metric(
			timepoints,
			detection_cutoff = options.detection_breakpoint,
			fixed_cutoff = options.fixed_breakpoint,
			metric = options.metric
		)
		PAIRWISE_CALCULATIONS.update(pair_array.copy())
	return PAIRWISE_CALCULATIONS


def generate_genotypes(timepoints: pandas.DataFrame, options: GenotypeOptions) -> Tuple[pandas.DataFrame, pandas.Series, Any]:
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

	pairwise_calculations = _update_pairwise_array(timepoints,options)
	genotypes, linkage_matrix = calculate_genotypes_from_given_method(
		timepoints,
		pairwise_calculations,
		options.method,
		options.similarity_breakpoint,
		options.difference_breakpoint,
		options.starting_genotypes
	)
	_mean_genotypes = calculate_mean_genotype(genotypes, timepoints)
	genotype_members = _mean_genotypes.pop('members')
	_mean_genotypes = _mean_genotypes[sorted(_mean_genotypes.columns)]

	return _mean_genotypes, genotype_members, linkage_matrix

def generate_genotypes_with_filter(original_timepoints:pandas.DataFrame, options:GenotypeOptions, frequency_breakpoints:List[float], strict_filter:bool):
	original_genotypes, original_genotype_members, linkage_matrix = generate_genotypes(original_timepoints, options)

	timepoints, mean_genotypes, genotype_members, linkage_matrix = filter_genotypes(
		original_timepoints,
		options,
		frequency_breakpoints,
		strict_filter
	)

	original_genotypes['members'] = genotype_members
	_tm = map_trajectories_to_genotype(original_genotype_members)
	original_timepoints['genotype'] = [_tm.get(i) for i in original_timepoints.index]
	return original_genotypes, timepoints, mean_genotypes, genotype_members, linkage_matrix



if __name__ == "__main__":
	pass
