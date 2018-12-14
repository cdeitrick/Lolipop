from typing import Dict, Tuple

import pandas

try:
	from muller.order_clusters import ClusterType
	from muller_genotypes import PairwiseArrayType
except ModuleNotFoundError:
	from order_clusters import ClusterType


def generate_ggmuller_population_table(mean_genotypes: pandas.DataFrame, edges: pandas.DataFrame, detection_cutoff: float) -> pandas.DataFrame:
	"""
		Converts the genotype frequencies to a population table suitable for ggmuller.
	Parameters
	----------
	mean_genotypes: pandas.DataFrame
		Contains the mean frequency for each genotype at each timepoint.
	edges: pandas.DataFrame
		The output from create_ggmuller_edges()
	detection_cutoff: float
		The cutoff to determine whether a trajectory or genotype counts as being nonzero at a given timepoint.
	Returns
	-------

	"""

	# "Generation", "Identity" and "Population"
	# Adjust populations to account for inheritance.
	# If a child genotype fixed, the parent genotype should be replaced.

	# Use a copy of the dataframe to avoid making changes to the original.
	modified_genotypes: pandas.DataFrame = mean_genotypes.copy(True)
	modified_genotypes = modified_genotypes[modified_genotypes.index.isin(edges['Identity'])]

	# Generate a list of all muller_genotypes that arise in the background of each genotype.
	children = dict()
	for _, row in edges.iterrows():
		parent = row['Parent']
		identity = row['Identity']

		children[parent] = children.get(parent, list()) + [identity]

	table = list()

	for genotype_label, genotype in modified_genotypes.iterrows():
		if genotype_label in children:
			genotype_frequencies = genotype[genotype > detection_cutoff]
			genotype_children: pandas.Series = modified_genotypes.loc[children[genotype_label]].max()
			genotype_frequencies: pandas.Series = genotype_frequencies - genotype_children
			genotype_frequencies = genotype_frequencies.mask(lambda s: s < detection_cutoff, 0.01)
			genotype_frequencies = genotype_frequencies.dropna()

		else:
			genotype_frequencies = genotype

		# Remove timepoints where the value was 0 or less than 0 due to above line.

		for timepoint, frequency in genotype_frequencies.items():
			row = {
				'Identity':   genotype_label,
				'Generation': timepoint,
				'Population': frequency * 100
			}
			table.append(row)

	temp_df = pandas.DataFrame(table)

	generation_groups = temp_df.groupby(by = 'Generation')
	for generation, group in generation_groups:
		population = group['Population'].sum()
		if population <= 100:
			p = 100 - population
		else:
			p = 0
		table.append({'Generation': generation, "Identity": "genotype-0", "Population": p})

	return pandas.DataFrame(table)


def generate_ggmuller_edges_table(genotype_clusters: ClusterType) -> pandas.DataFrame:
	table = list()
	for genotype_info in genotype_clusters.values():
		identity = genotype_info.name
		background = genotype_info.background

		if len(background) == 1:
			parent = 'genotype-0'
		else:
			parent = background[0]
		if parent == identity:
			parent = 'genotype-0'

		row = {
			'Parent':   parent,
			'Identity': identity
		}
		table.append(row)
	table = pandas.DataFrame(table)[['Parent', 'Identity']]

	return table


def generate_p_value_table(p_values, trajectory_genotypes: Dict[str, str]) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	timeseries_table = list()
	for (left, right), calculation in p_values.items():
		if calculation.mean_series is None:
			# Assume all tables are None
			continue

		pair_df = pandas.DataFrame([calculation.sigma_series, calculation.mean_series, calculation.difference_series])
		pair_df['variable'] = ['sigma', 'mean', 'difference']
		pair_df['leftTrajectory'] = left
		pair_df['rightTrajectory'] = right
		pair_df['sigmaPair'] = calculation.sigma
		pair_df['differencePair'] = calculation.difference_mean
		pair_df['pvalue'] = calculation.pvalue

		timeseries_table.append(pair_df)

	df = pandas.concat(timeseries_table, sort = True)
	numeric_columns = [i for i in df.columns if isinstance(i, int)]
	nonnumeric_columns = [i for i in df.columns if i not in numeric_columns]
	df = df[nonnumeric_columns + numeric_columns]
	matrix = generate_p_value_matrix(p_values, trajectory_genotypes)

	return df, matrix


def generate_p_value_matrix(p_values, trajectory_genotypes: Dict[str, str]):
	""" Converts a dictionary mapping pairs of trajectories with thier respoective p-values into a similarity matrix."""
	import itertools
	import math

	p_values = {k: v.pvalue for k, v in p_values.items()}
	table = list()
	all_ids = sorted(
		set(itertools.chain.from_iterable(p_values.keys())),
		key = lambda s: (trajectory_genotypes.get(s, 'zzz'), s)
	)

	for index in all_ids:
		row = dict()
		for column in all_ids:
			if column == index:
				value = math.nan
			else:
				value = p_values[column, index]
			row[column] = value
		series = pandas.Series(row, name = index)
		table.append(series)
	df = pandas.DataFrame(table)
	return df


def generate_missing_trajectories_table(trajectories: pandas.DataFrame, original_trajectories: pandas.DataFrame) -> pandas.DataFrame:
	missing_trajectories = original_trajectories[~original_trajectories.index.isin(trajectories.index)]
	concat_trajectories = pandas.concat([trajectories, missing_trajectories], sort = False)

	return concat_trajectories


def generate_trajectory_table(trajectories: pandas.DataFrame, parent_genotypes: pandas.Series, info: pandas.DataFrame) -> pandas.DataFrame:
	# Sorts the trajectory table and adds additional columns from the original table.
	trajectories = trajectories[sorted(trajectories.columns, key = lambda s: int(s))]
	if info is not None:
		trajectory_table: pandas.DataFrame = trajectories.copy()
		trajectory_table['genotype'] = [parent_genotypes[k] for k in trajectory_table.index]
		trajectory_table = trajectory_table.join(info).sort_values(by = ['genotype'])
		return trajectory_table


if __name__ == "__main__":
	pass
