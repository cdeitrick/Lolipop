from typing import List, Optional, Tuple

import numpy
import pandas
from loguru import logger

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


def _update_pairwise_array(timepoints: pandas.DataFrame, dlimit:float, flimit:float, metric:str):
	global PAIRWISE_CALCULATIONS
	logger.debug(f"Pairwise calculations already exist: {bool(PAIRWISE_CALCULATIONS)}")
	if PAIRWISE_CALCULATIONS:
		# PAIRWISE_CALCULATIONS already contains previous calcuations.
		_before = len(PAIRWISE_CALCULATIONS)
		PAIRWISE_CALCULATIONS = PAIRWISE_CALCULATIONS.reduce(timepoints.index)
		_after = len(PAIRWISE_CALCULATIONS)
		logger.debug(f"Reduced Pairwise calculations from {_before} to {_after}")
	else:
		pair_array = calculate_pairwise_metric(
			timepoints,
			detection_cutoff = dlimit,
			fixed_cutoff = flimit,
			metric = metric
		)
		PAIRWISE_CALCULATIONS.update(pair_array.copy())
	return PAIRWISE_CALCULATIONS


def calculate_genotypes(timepoints: pandas.DataFrame, dlimit:float, flimit:float, sbreakpoint:float, dbreakpoint:float, method:str, metric:str, starting_genotypes:List[List[str]]) -> Tuple[pandas.DataFrame, pandas.Series, Optional[numpy.array]]:
	pairwise_calculations = _update_pairwise_array(timepoints, dlimit, flimit, metric)
	genotypes, linkage_matrix = calculate_genotypes_from_given_method(
		timepoints,
		pairwise_calculations,
		method,
		sbreakpoint,
		dbreakpoint,
		starting_genotypes
	)
	mean_genotypes = calculate_mean_genotype(genotypes, timepoints)
	genotype_members = mean_genotypes.pop('members')
	return mean_genotypes, genotype_members, linkage_matrix


def generate_genotypes(trajectories: pandas.DataFrame, dlimit: float, flimit: float, similarity_breakpoint:float, difference_breakpoint: float, method: str,
		metric: str, breakpoints: List[float] = None, starting_genotypes: List[List[str]] = None) -> Tuple[
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
	dlimit, flimit:float
		The detection, significant, and fixed cutoff values, respectively
	similarity_breakpoint: float
		The cutoff value to use when clustering trajectories into genotypes. Two trajectories must have a distance less than this value to be
		considered members of the same genotype.
	difference_breakpoint: float
		Only used when using the two-step method. Governs when the genotype contains two trajectories which should be split into separate genotypes.
	method: {'two-step', 'hierarchy'}
		The clustering method to use. The hierarchical method uses the maximum distance between a point and a candidate cluster to determine whether
		to group said point to the cluster.
	metric: {'binomial', 'pearson', 'minkowski'}
		The distance metric to determine trajectory similarity.
	breakpoints: Optional[List[float]]
		An empty/None breakpoints value indicated that the filters should be skipped.
	starting_genotypes: List[List[str]]
		A list of genotypes to start the clustering algorithm with. The distance metrics will be modified so that the trajectories specified are
		grouped together.
	Returns
	-------
	pandas.DataFrame, pandas.Series, numpy.array
		- The genotype table
		- A map of genotypes to members.
		- A linkage matrix, if hierarchical clustering was used.
	"""
	use_strict_filter = False
	modified_trajectories = trajectories.copy(deep = True)  # To avoid unintended changes
	genotype_table, genotype_members, linkage_matrix = calculate_genotypes(
		modified_trajectories,
		dlimit, flimit,
		similarity_breakpoint, difference_breakpoint,
		method, metric,
		starting_genotypes
	)

	if breakpoints:
		_iterations = 20  # arbitrary, used to make sure the program does not encounter an infinite loop.
	else:
		_iterations = 0  # The for loop shouldn't excecute.

	for index in range(_iterations):
		invalid_members = filters.filter_genotypes(genotype_table, genotype_members, dlimit, breakpoints, use_strict_filter)
		if invalid_members:
			# Remove these trajectories from the trajectories table.
			modified_trajectories = modified_trajectories[~modified_trajectories.index.isin(invalid_members)]
			# Re-calculate the genotypes based on the remaining trajectories.
			genotype_table, genotype_members, linkage_matrix = calculate_genotypes(modified_trajectories,
				dlimit, flimit,
			similarity_breakpoint, difference_breakpoint,
			method, metric,
			starting_genotypes
			)
		else:
			break
	return genotype_table, genotype_members, linkage_matrix


if __name__ == "__main__":
	pass
