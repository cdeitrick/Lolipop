import logging
from typing import List, Optional, Tuple

import numpy
import pandas

from options import GenotypeOptions

logger = logging.getLogger(__name__)
try:
	from muller.clustering.average import calculate_mean_genotype
	from clustering.metrics import PairwiseCalculationCache, calculate_pairwise_metric
	from clustering.methods import calculate_genotypes_from_given_method
	from clustering import filters
except ModuleNotFoundError:
	from .average import calculate_mean_genotype
	from .metrics import PairwiseCalculationCache, calculate_pairwise_metric
	from .methods import calculate_genotypes_from_given_method
	from . import filters

PAIRWISE_CALCULATIONS = PairwiseCalculationCache()


def _update_pairwise_array(timepoints: pandas.DataFrame, options: GenotypeOptions):
	global PAIRWISE_CALCULATIONS
	logger.info(f"Pairwise calculations already exist: {bool(PAIRWISE_CALCULATIONS)}")
	if PAIRWISE_CALCULATIONS:
		# PAIRWISE_CALCULATIONS already contains previous calcuations.
		_before = len(PAIRWISE_CALCULATIONS)
		PAIRWISE_CALCULATIONS = PAIRWISE_CALCULATIONS.reduce(timepoints.index)
		_after = len(PAIRWISE_CALCULATIONS)
		logging.info(f"Reduced Pairwise calculations from {_before} to {_after}")
	else:
		pair_array = calculate_pairwise_metric(
			timepoints,
			detection_cutoff = options.detection_breakpoint,
			fixed_cutoff = options.fixed_breakpoint,
			metric = options.metric
		)
		PAIRWISE_CALCULATIONS.update(pair_array.copy())
	return PAIRWISE_CALCULATIONS


def calculate_genotypes(timepoints: pandas.DataFrame, options: GenotypeOptions) -> Tuple[pandas.DataFrame, pandas.Series, Optional[numpy.array]]:
	pairwise_calculations = _update_pairwise_array(timepoints, options)
	genotypes, linkage_matrix = calculate_genotypes_from_given_method(
		timepoints,
		pairwise_calculations,
		options.method,
		options.similarity_breakpoint,
		options.difference_breakpoint,
		options.starting_genotypes
	)
	mean_genotypes = calculate_mean_genotype(genotypes, timepoints)
	genotype_members = mean_genotypes.pop('members')
	return mean_genotypes, genotype_members, linkage_matrix


def generate_genotypes(trajectories: pandas.DataFrame, options: GenotypeOptions, breakpoints: List[float] = None) -> Tuple[
	pandas.DataFrame, pandas.Series, Optional[numpy.array]]:
	"""
	Parameters
	----------
	trajectories: pandas.DataFrame
		A timeseries dataframe, usually generated from `import_table.import_trajectory_table`.
			- Index -> str
				Names unique to each trajectory.
			- Columns -> int
				The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	options: GenotypeOptions
		An instance of `GenotypeOptions` with the desired options. Can be generated automatically via `GenotypeOptions.from_breakpoints(0.03)`.
	breakpoints: Optional[List[float]]
		An empty/None breakpoints value indicated that the filters should be skipped.
	Returns
	-------
	pandas.DataFrame, pandas.Series, numpy.array
		- The genotype table
		- A map of genotypes to members.
		- A linkage matrix, if hierarchical clustering was used.
	"""
	use_strict_filter = False
	modified_trajectories = trajectories.copy(deep = True)  # To avoid unintended changes
	genotype_table, genotype_members, linkage_matrix = calculate_genotypes(modified_trajectories, options)

	if breakpoints:
		_iterations = 20  # arbitrary, used to make sure the program does not encounter an infinite loop.
	else:
		_iterations = 0  # The for loop shouldn't excecute.

	for index in range(_iterations):
		invalid_members = filters.filter_genotypes(genotype_table, genotype_members, options, breakpoints, use_strict_filter)
		if invalid_members:
			# Remove these trajectories from the trajectories table.
			modified_trajectories = modified_trajectories[~modified_trajectories.index.isin(invalid_members)]
			# Re-calculate the genotypes based on the remaining trajectories.
			genotype_table, genotype_members, linkage_matrix = calculate_genotypes(modified_trajectories, options = options)
		else:
			break
	return genotype_table, genotype_members, linkage_matrix


if __name__ == "__main__":
	pass
