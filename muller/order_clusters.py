from pathlib import Path
import pandas
import math
import itertools

DETECTION_THRESHOLD = 0.03
FREQUENCY_THRESHOLD = 0.15
FIXED_THRESHOLD = 0.85

FREQUENCY_BREAKPOINTS = [0.90, 0.75, 0.60, 0.45, 0.30, 0.15, 0.00]  # Used when sorting non-fixed genotypes.


def sort_genotypes(genotype_trajectories: pandas.DataFrame) -> pandas.DataFrame:
	# sorted_df = sort_genotype_frequencies(genotype_trajectories, FIXED_THRESHOLD)
	sorted_genotypes = list()
	genotypes: pandas.DataFrame = genotype_trajectories.copy()
	for frequency in [FIXED_THRESHOLD] + FREQUENCY_BREAKPOINTS:
		sorted_df = sort_genotype_frequencies(genotypes, frequency)
		if not sorted_df.empty:
			genotypes = genotypes.drop(sorted_df.index)
			sorted_genotypes.append(sorted_df)

	df = pandas.concat(sorted_genotypes)
	df = genotype_trajectories.reindex(df.index)
	return df


def sort_genotype_frequencies(genotype_trajectories: pandas.DataFrame, frequency_breakpoint) -> pandas.DataFrame:
	# Should be 1, 7, 2|3, 6, 8|4, 5, 14, 9|17, 19|13|..., 21|10|15, 18|11
	# Should Be 8, 7|12|28,
	#frequency_threshold = FIXED_THRESHOLD if FREQUENCY_THRESHOLD < frequency_breakpoint else frequency_breakpoint
	sorted_genotypes = list()
	transposed = genotype_trajectories.transpose()

	first_detected: pandas.Series = (transposed > DETECTION_THRESHOLD).idxmax(0).sort_values()
	first_detected.name = 'firstDetected'

	first_above_threshold: pandas.Series = (transposed > FREQUENCY_THRESHOLD ).idxmax(0).sort_values()
	first_above_threshold.name = 'firstThreshold'

	first_fixed: pandas.Series = (transposed > frequency_breakpoint).idxmax(0).sort_values()
	first_fixed.name = 'firstFixed'

	# Remove the genotypes which were never detected or never rose above the threshold.
	first_detected_reduced = first_detected.iloc[first_detected.nonzero()]
	first_above_threshold_reduced = first_above_threshold.replace(0, 13) # To replicate the behavior in the matlab script
	first_fixed_reduced = first_fixed.iloc[first_fixed.nonzero()]

	df: pandas.DataFrame = pandas.concat([first_fixed_reduced, first_detected_reduced, first_above_threshold_reduced],
		axis = 1)
	df = df[~df['firstFixed'].isna()]

	if frequency_breakpoint == FIXED_THRESHOLD:
		df = df.dropna()
	sorted_df = df.sort_values(by = ['firstDetected', 'firstThreshold'])

	# Sort genotypes based on frequency if two or more share the same fixed timepoint, detection timpoint, threshold timepoint.

	debug = False
	if not debug:
		freq_groups = list()
		groups = sorted_df.groupby(by = list(sorted_df.columns))
		for label, group in groups:
			trajectories = genotype_trajectories.loc[group.index]
			if len(group) < 2:
				freq_groups.append(trajectories)
			else:
				trajectories = genotype_trajectories.drop([i for i in genotype_trajectories.index if i not in group.index])
				#sgens = gens.sort_values(by = list(genotype_trajectories.columns))
				freq_groups.append(trajectories)
		if freq_groups: # Make sure freq_groups is not empty
			freq_df = pandas.concat(freq_groups)
		else:
			freq_df = sorted_df
	else:
		freq_df = sorted_df
	return freq_df


def order_clusters(genotype_trajectories: pandas.DataFrame):
	"""

	Parameters
	----------
	genotype_trajectories: pandas.DataFrame
		Output from get_mean_genotypes

	Returns
	-------

	"""
	sort_genotypes(genotype_trajectories)
	for timepoint, trajectories in genotype_trajectories.iteritems():
		fixed_trajectories: pandas.Series = trajectories[trajectories > 0.85]  # if there are fixed trajectories at time

		# arrange all fixed genotypes by when they were first detected and first above 15%
		if not fixed_trajectories.empty:
			sorted_trajectories = fixed_trajectories.sort_values()


if __name__ == "__main__":
	from muller import variables
	from muller.time_series_import import import_timeseries
	from muller.get_genotypes import get_genotypes, get_mean_genotypes

	filename = Path(variables.filename)
	timepoints, info = import_timeseries(filename)

	genotypes = get_genotypes(timepoints)
	mean_genotypes = get_mean_genotypes(genotypes, timepoints)

	order_clusters(mean_genotypes)
