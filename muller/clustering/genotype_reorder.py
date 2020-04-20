import pandas

try:
	from muller.inheritance import timepoint_detection
except ModuleNotFoundError:
	from muller.inheritance import timepoint_detection


class SortGenotypeTableWorkflow:
	"""
		Sorts the muller_genotypes based on when they were first detected and first fixed.
		dlimit, slimit, flimit: float
			The three breakpoints that are used to determine the sort order of the genotypes.
		breakpoints: List[float]
			Frequencies with which to group genotypes into sortable bins. Each bin will be sorted individually then added to the final output.
	"""

	def __init__(self, dlimit: float, flimit: float):
		self.dlimit = dlimit
		self.slimit = 0.15
		self.flimit = flimit
		# The `breakpoints` value is a bit arbitrary, so it should be safe to hard-code it.
		# 	This will actually prevent the most common error when sorting genotypes (i.e. no breakpoints given) so it's worth
		#	hard-coding it to prevent that issue.
		# TODO: Make sure this is hard-coded globally.
		# self.breakpoints = sorted([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0] + [self.dlimit, self.flimit], reverse = True)
		self.breakpoints = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

	@staticmethod
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
			# if ft == 130:  # dummy value assigned above. Revert back to 0 since '130' isn't a valid timepoint.
			#	ft = 0
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

	def _sort_genotype_frequencies(self, genotype_trajectories: pandas.DataFrame, frequency_breakpoint: float) -> pandas.DataFrame:
		"""
			Sorts the genotype table based on the timepoints that first satisfy the threshold requirements.
		Parameters
		----------
		genotype_trajectories
		frequency_breakpoint
		"""
		# Get the values for each genotype at the first timepoint. This will be used to check if the genotype existed at the first timepoint
		# or has an initial detected timepoint of 0 due to the way pandas returns indicies.
		detection_cutoff = min(self.dlimit, frequency_breakpoint)
		transposed = genotype_trajectories.transpose()

		first_detected_reduced = timepoint_detection.get_first_detected_timepoint(transposed, detection_cutoff)
		first_above_threshold_reduced = timepoint_detection.get_first_significant_timepoint(transposed, self.slimit)
		# Use the frequency breakpoint rather than the fixed cutoff.
		first_fixed_reduced = timepoint_detection.get_first_fixed_timepoint(transposed, frequency_breakpoint)

		combined_df: pandas.DataFrame = pandas.concat(
			[first_fixed_reduced, first_detected_reduced, first_above_threshold_reduced],
			axis = 1, sort = False)
		df = combined_df[~combined_df['firstFixed'].isna()]

		if frequency_breakpoint == self.flimit:
			df = df.dropna()
		df = df[['firstDetected', 'firstFixed', 'firstThreshold']]
		sorted_frequencies = df.sort_values(by = ['firstDetected', 'firstThreshold'], ascending = True)

		# Sort genotypes based on frequency if two or more share the same fixed timepoint, detection timpoint, threshold timepoint.
		# Iterate over the conbinations of 'firstFixed', 'firstDetected', and 'firstThreshold' and sort trajectories that belong to the sample combination.
		# Need to build a new trajectory table based on the sorted timepoint table.
		freq_df = self._build_sorted_frequency_table(genotype_trajectories, sorted_frequencies)

		return freq_df

	@staticmethod
	def remove_trajectories_below_threshold(df: pandas.DataFrame, threshold: float) -> pandas.DataFrame:
		""" Removes all genotypes that do not have at least one frequency above the given threshold."""
		result = df[df.max(axis = 1) >= threshold]
		return result

	def run(self, unsorted_genotypes: pandas.DataFrame):
		"""
			Sorts the muller_genotypes based on when they were first detected and first fixed.
		Parameters
		----------
		unsorted_genotypes: pandas.Dataframe
			A dataframe with the mean frequency of each genotype, derived from the member trajectories
			in that genotype. Each row should correspond to a single genotype.
		Returns
		-------
		sorted_genotype: pandas.DataFrame
		"""
		sorted_genotypes = list()
		current_genotypes: pandas.DataFrame = unsorted_genotypes.copy()
		for frequency in self.breakpoints:
			# Ignore genotypes that do not have at least on timepoint exceeding the current frequency.
			# This should be taken care of suring the filtering step.
			genotypes_above_threshold = self.remove_trajectories_below_threshold(current_genotypes, frequency)

			sorted_dataframe = self._sort_genotype_frequencies(
				genotype_trajectories = genotypes_above_threshold,
				frequency_breakpoint = frequency
			)
			if sorted_dataframe is not None:
				current_genotypes = current_genotypes.drop(sorted_dataframe.index)
				sorted_genotypes.append(sorted_dataframe)
		df = pandas.concat(sorted_genotypes, sort = False)
		df = unsorted_genotypes.reindex(df.index)

		# Make sure the genotype column is labelled correctly.
		df.index.name = "Genotype"

		return df
