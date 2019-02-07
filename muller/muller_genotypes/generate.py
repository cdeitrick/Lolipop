import logging
from typing import Any, List, Tuple

import pandas

from options import GenotypeOptions

logger = logging.getLogger(__name__)
try:
	from muller.muller_genotypes.average import calculate_mean_genotype
	from muller_genotypes.metrics import PairwiseCalculation, calculate_pairwise_metric
	from muller_genotypes.methods import calculate_genotypes_from_given_method
	from muller_genotypes.filters import get_fuzzy_backgrounds, filter_trajectories, find_first_invalid_genotype
	from muller.widgets import map_trajectories_to_genotype
except ModuleNotFoundError:
	from .average import calculate_mean_genotype
	from .metrics import PairwiseCalculation, calculate_pairwise_metric
	from .methods import calculate_genotypes_from_given_method
	from .filters import get_fuzzy_backgrounds, filter_trajectories, find_first_invalid_genotype
	from widgets import map_trajectories_to_genotype

PAIRWISE_CALCULATIONS = PairwiseCalculation()


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

	pairwise_calculations = _update_pairwise_array(timepoints, options)
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


def generate_genotypes_with_filter(original_timepoints: pandas.DataFrame, options: GenotypeOptions, frequency_breakpoints: List[float],
		strict_filter: bool):
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


def filter_genotypes(trajectory_table: pandas.DataFrame, goptions: GenotypeOptions, frequency_cutoffs: List[float],
		use_strict_filter: bool) -> Tuple[pandas.DataFrame, pandas.DataFrame, Any, Any]:
	"""
		Iteratively calculates the population genotypes, checks and removes invalid genotypes, and recomputes the genotypes until no changes occur.
	Parameters
	----------
	trajectory_table: pandas.DataFrame
	goptions: GenotypeOptions
	frequency_cutoffs: List[float]
	use_strict_filter: bool

	Returns
	-------

	"""
	# Remove 1.0 fro mthe list of frequency breakpoints to account for measurement errors.
	logger.info("Filtering genotypes...")
	frequency_cutoffs = [i for i in frequency_cutoffs if i <= goptions.fixed_breakpoint]
	logger.info("\tFrequency cutoffs:" + str(frequency_cutoffs))
	trajectory_table = trajectory_table.copy(deep = True)  # To avoid any unintended changes to the original table.
	filtered_trajectory_table = filter_trajectories(trajectory_table, goptions.detection_breakpoint, goptions.fixed_breakpoint)

	# Generate the initial genotypes.
	genotype_table, genotype_members, linkage_table = generate_genotypes(filtered_trajectory_table, options = goptions)

	_iterations = 20  # arbitrary, used to ensure the program does not encounter an infinite loop.
	for index in range(_iterations):
		logger.info(f"filtering iteration {index} of {_iterations}")
		# Find all the backgrounds for this population. Some may fall below the usual `fixed_cutoff` threshold, so use the same frequency breakpoints
		# used when sorting the genotypes.
		current_backgrounds, (dlimit, flimit) = get_fuzzy_backgrounds(genotype_table, frequency_cutoffs)
		logger.info(f"Backgrounds:" + str(list(current_backgrounds.index)))

		# Search for genotypes that do not make sense in the context of an evolved population.
		fuzzy_detected_cutoff = max(goptions.detection_breakpoint, dlimit)
		logger.info(f"fuzzy cutoffs: {fuzzy_detected_cutoff}, {flimit}")
		current_invalid_genotype = find_first_invalid_genotype(genotype_table, current_backgrounds, fuzzy_detected_cutoff, flimit, use_strict_filter)
		logger.info(f"Invalid genotype: {current_invalid_genotype}")
		if current_invalid_genotype is None:
			break
		else:
			# Get a list of the trajectories that form this genotype.
			invalid_members = genotype_members.loc[current_invalid_genotype].split('|')
			logger.info("Invalid members: " + str(invalid_members))
			# Remove these trajectories from the trajectories table.
			filtered_trajectory_table = filtered_trajectory_table[~filtered_trajectory_table.index.isin(invalid_members)]
			# Re-calculate the genotypes based on the remaining trajectories.
			genotype_table, genotype_members, linkage_table = generate_genotypes(filtered_trajectory_table, options = goptions)

	# cache.append((trajectory_table.copy(), genotype_table.copy()))
	# Update the trajectories that comprise each genotype.
	else:
		print(f"Could not filter the genotypes after {_iterations} iterations.")
	return filtered_trajectory_table, genotype_table, genotype_members, linkage_table


if __name__ == "__main__":
	pass
