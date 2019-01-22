import pandas

from options import SortOptions


def sort_genotypes(genotype_frequencies: pandas.DataFrame, options: SortOptions) -> pandas.DataFrame:
	"""
		Sorts the muller_genotypes based on when they were first detected and first fixed.
	Parameters
	----------
	genotype_frequencies:pandas.Dataframe
		A dataframe with the mean frequency of each genotype, derived from the member trajectories
		in that genotype. Each row should correspond to a single genotype.
	options: SortOptions

	Returns
	-------
	sorted_genotype: pandas.DataFrame
	"""
	sorted_genotypes = list()
	current_genotypes: pandas.DataFrame = genotype_frequencies.copy()
	for frequency in [options.fixed_breakpoint] + options.frequency_breakpoints:
		sorted_dataframe = sort_genotype_frequencies(
			genotype_trajectories = current_genotypes,
			frequency_breakpoint = frequency,
			detection_cutoff = options.detection_breakpoint,
			significant_cutoff = options.significant_breakpoint,
			fixed_cutoff = options.fixed_breakpoint
		)

		if not sorted_dataframe.empty:
			current_genotypes = current_genotypes.drop(sorted_dataframe.index)
			sorted_genotypes.append(sorted_dataframe)

	df = pandas.concat(sorted_genotypes, sort = False)
	df = genotype_frequencies.reindex(df.index)
	return df


def sort_genotype_frequencies(genotype_trajectories: pandas.DataFrame, frequency_breakpoint: float,
		detection_cutoff: float, significant_cutoff: float, fixed_cutoff: float) -> pandas.DataFrame:
	transposed = genotype_trajectories.transpose()

	# Finds the first timepoint where each genotype exceeded the fixed frequency
	first_detected: pandas.Series = (transposed > detection_cutoff).idxmax(0).sort_values()
	first_detected.name = 'firstDetected'

	first_above_threshold: pandas.Series = (transposed > significant_cutoff).idxmax(0).sort_values()
	first_above_threshold.name = 'firstThreshold'

	# noinspection PyUnresolvedReferences
	first_fixed: pandas.Series = (transposed > frequency_breakpoint).idxmax(0).sort_values()
	first_fixed.name = 'firstFixed'

	# Remove the muller_genotypes which were never detected or never rose above the threshold.
	first_detected_reduced = first_detected.iloc[first_detected.nonzero()]
	# To replicate the behavior in the matlab script

	first_above_threshold_reduced = first_above_threshold.replace(0, 130)

	first_fixed_reduced: pandas.DataFrame = first_fixed.iloc[first_fixed.nonzero()]

	# if first_fixed_reduced.empty:
	#	first_fixed_reduced = first_fixed

	combined_df: pandas.DataFrame = pandas.concat(
		[first_fixed_reduced, first_detected_reduced, first_above_threshold_reduced],
		axis = 1, sort = False)

	df = combined_df[~combined_df['firstFixed'].isna()]

	if frequency_breakpoint == fixed_cutoff:
		df = df.dropna()
	sorted_frequencies = df.sort_values(by = ['firstDetected', 'firstThreshold'])

	# Sort genotypes based on frequency if two or more share the same fixed timepoint, detection timpoint, threshold timepoint.
	freq_groups = list()
	groups = sorted_frequencies.groupby(by = list(sorted_frequencies.columns))

	# Iterate over the conbinations of 'firstFixed', 'firstDetected', and 'firstThreshold' and sort trajectories that belong to the sample combination.
	for (ff, fd, ft), group in groups:
		if ft == 130:  # dummy value assigned above. Revert back to 0 since '130' isn't a valid timepoint.
			ft = 0
		trajectories: pandas.DataFrame = genotype_trajectories.loc[group.index]

		if len(group) < 2:
			# There is only one genotype in this combination of timepoints. No need to sort.
			freq_groups.append(trajectories)
		else:
			# More than one genotype share this combination of key timepoints. Sort by frequency.
			# Sort from highest to lowest
			trajectories = trajectories.sort_values(by = [ff, fd, ft], ascending = False)
			freq_groups.append(trajectories)

	if freq_groups:  # Make sure freq_groups is not empty
		freq_df = pandas.concat(freq_groups)
	else:
		freq_df = sorted_frequencies

	return freq_df


if __name__ == "__main__":
	pass
