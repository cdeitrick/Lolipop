from typing import List

import pandas


def _calculate_mean_frequencies_of_trajectories(name: str, genotype_timeseries: pandas.DataFrame, genotype: List[str]) -> pandas.Series:
	mean_genotype_timeseries = genotype_timeseries.mean()
	mean_genotype_timeseries['members'] = "|".join(map(str, genotype))
	mean_genotype_timeseries.name = name

	return mean_genotype_timeseries


def calculate_mean_genotype(all_genotypes: List[List[str]], timeseries: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Calculates the mean frequency of each genotype ate every timepoint.
	Parameters
	----------
	all_genotypes: List[List[str]
		A list of all muller_genotypes for a given population
	timeseries: pandas.DataFrame
		Must have a 'Trajectory' column along with the columns of the original table that represent timepoints.
		Each row corresponds to a single mutational trajectory.

	Returns
	-------
	A dataframe where every row corresponds to a genotype.
	member trajectories are listed under the 'members' column.
	every column represents a timepoint.
	"""
	mean_genotypes = list()
	for index, genotype in enumerate(all_genotypes, start = 1):
		genotype_timeseries = timeseries.loc[genotype]
		mean_genotype_timeseries = _calculate_mean_frequencies_of_trajectories(f"genotype-{index}", genotype_timeseries, genotype)
		mean_genotypes.append(mean_genotype_timeseries)

	mean_genotypes = pandas.DataFrame(mean_genotypes)
	# For consistency
	mean_genotypes.index.name = 'Genotype'

	return mean_genotypes
