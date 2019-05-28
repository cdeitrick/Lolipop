from typing import List

import pandas
from loguru import logger

try:
	from muller.inheritance import timepoint_detection
except ModuleNotFoundError:
	from . import timepoint_detection


def sort_genotypes(genotype_frequencies: pandas.DataFrame, dlimit: float, slimit: float, flimit: float, breakpoints: List[float]) -> pandas.DataFrame:
	"""
		Sorts the muller_genotypes based on when they were first detected and first fixed.
	Parameters
	----------
	genotype_frequencies:pandas.Dataframe
		A dataframe with the mean frequency of each genotype, derived from the member trajectories
		in that genotype. Each row should correspond to a single genotype.
	dlimit, slimit, flimit: float
		The three breakpoints that are used to determine the sort order of the genotypes.
	breakpoints: List[float]
		Frequencies with which to group genotypes into sortable bins. Each bin will be sorted individually then added to the final output.

	Returns
	-------
	sorted_genotype: pandas.DataFrame
	"""
	sorted_genotypes = list()
	current_genotypes: pandas.DataFrame = genotype_frequencies.copy()
	for frequency in [flimit] + breakpoints:
		logger.debug(f"filtering based on frequency {frequency}")
		# Ignore genotypes that do not have at least on timepoint exceeding the current frequency.
		genotypes_above_threshold = _remove_low_frequency_series(current_genotypes, frequency)

		sorted_dataframe = _sort_genotype_frequencies(
			genotype_trajectories = genotypes_above_threshold,
			frequency_breakpoint = frequency,
			detection_cutoff = dlimit,
			significant_cutoff = slimit,
			fixed_cutoff = flimit
		)
		if sorted_dataframe is not None:
			current_genotypes = current_genotypes.drop(sorted_dataframe.index)
			sorted_genotypes.append(sorted_dataframe)
	df = pandas.concat(sorted_genotypes, sort = False)
	df = genotype_frequencies.reindex(df.index)
	return df


def _remove_low_frequency_series(df: pandas.DataFrame, threshold: float) -> pandas.DataFrame:
	""" Removes all genotypes that do not have at least one frequency above the given threshold."""
	result = df[df.max(axis = 1) >= threshold]
	return result


# noinspection PyUnresolvedReferences,PyTypeChecker
def _get_timepoint_above_threshold(transposed_timepoints: pandas.DataFrame, cutoff: float, name: str = None) -> pandas.Series:
	"""
		Calculates when a genotype was first fixed based on `cutoff`. The resulting series will be named `name`, if given.
	Parameters
	----------
	transposed_timepoints: pandas.DataFrame
		A dataframe indexed by timepoint and with each series a separate column.
	cutoff: float
		The threshold a frequency must exceed to be considered a valid timepoint.

	Returns
	-------
	pandas.Series
		A series mapping series name to the first timepoint that series first exceeded the threshold value.
		Series that never exceed the threshold are mapped to a value of 0.
	"""
	threshold_series: pandas.Series = (transposed_timepoints > cutoff).idxmax(0).sort_values()
	if name:
		threshold_series.name = name
	return threshold_series


def _sort_genotype_frequencies(genotype_trajectories: pandas.DataFrame, frequency_breakpoint: float,
		detection_cutoff: float, significant_cutoff: float, fixed_cutoff: float) -> pandas.DataFrame:
	"""
		Sorts the genotype table based on the timepoints that first satisfy the threshold requirements.
	Parameters
	----------
	genotype_trajectories
	frequency_breakpoint
	detection_cutoff
	significant_cutoff
	fixed_cutoff

	Returns
	-------
	"""
	# Get the values for each genotype at the first timepoint. This will be used to check if the genotype existed at the first timepoint
	# or has an initial detected timepoint of 0 due to the way pandas returns indicies.
	detection_cutoff = min(detection_cutoff, frequency_breakpoint)
	transposed = genotype_trajectories.transpose()

	first_detected_reduced = timepoint_detection.get_first_detected_timepoint(transposed, detection_cutoff)
	first_above_threshold_reduced = timepoint_detection.get_first_significant_timepoint(transposed, significant_cutoff)
	# Use the frequency breakpoint rather than the fixed cutoff.
	first_fixed_reduced = timepoint_detection.get_first_fixed_timepoint(transposed, frequency_breakpoint)

	combined_df: pandas.DataFrame = pandas.concat(
		[first_fixed_reduced, first_detected_reduced, first_above_threshold_reduced],
		axis = 1, sort = False)
	df = combined_df[~combined_df['firstFixed'].isna()]

	if frequency_breakpoint == fixed_cutoff:
		df = df.dropna()
	sorted_frequencies = df.sort_values(by = ['firstDetected', 'firstThreshold'], ascending = False)

	# Sort genotypes based on frequency if two or more share the same fixed timepoint, detection timpoint, threshold timepoint.
	# Iterate over the conbinations of 'firstFixed', 'firstDetected', and 'firstThreshold' and sort trajectories that belong to the sample combination.
	# Need to build a new trajectory table based on the sorted timepoint table.
	freq_df = _build_sorted_frequency_table(genotype_trajectories, sorted_frequencies)

	return freq_df


def _build_sorted_frequency_table(original_frequencies: pandas.DataFrame, thresholds: pandas.DataFrame) -> pandas.DataFrame:
	# Sort genotypes based on frequency if two or more share the same fixed timepoint, detection timpoint, threshold timepoint.
	# Iterate over the conbinations of 'firstFixed', 'firstDetected', and 'firstThreshold' and sort trajectories that belong to the sample combination.
	# Need to build a new trajectory table based on the sorted timepoint table.

	freq_groups = list()
	# Group the
	groups = thresholds.groupby(by = list(thresholds.columns))
	for (ff, fd, ft), group in groups:
		# Each group will look like this:
		#   	firstFixed firstDetected firstThreshold
		# 15          0            25              0
		# 10          0            25              0
		if ft == 130:  # dummy value assigned above. Revert back to 0 since '130' isn't a valid timepoint.
			ft = 0
		# Get a table of the original frequencies at each timepoint for the genotypes present in `group`
		trajectories: pandas.DataFrame = original_frequencies.loc[group.index]
		if len(group) < 2:
			# There is only one genotype in this combination of timepoints. No need to sort.
			freq_groups.append(trajectories)
		else:
			# More than one genotype share this combination of key timepoints. Sort by frequency.
			# Sort from highest to lowest using the timpoint columns as the sorting keys.
			# trajectories = trajectories.sort_values(by = [ff, ft, fd], ascending = False)
			trajectories = trajectories.sort_values(by = [fd, ft, ff], ascending = False)
			freq_groups.append(trajectories)

	try:
		freq_df = pandas.concat(freq_groups)
	except ValueError:
		# Can't concatenate an empty dataframe
		freq_df = None

	return freq_df
