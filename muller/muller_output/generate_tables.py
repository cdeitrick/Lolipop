from typing import Dict, List, Tuple, Union

import pandas


def _compile_parent_linkage(edges: pandas.DataFrame) -> Dict[str, List[str]]:
	""" Maps a genotype to a list of all genotypes that inherit from it"""
	children = dict()
	for _, row in edges.iterrows():
		parent = row['Parent']
		identity = row['Identity']

		children[parent] = children.get(parent, list()) + [identity]
	return children


def _subtract_children_from_parent(modified_genotypes: pandas.DataFrame, children: Dict[str, List[str]], detection_cutoff: float) -> pandas.DataFrame:
	children_table = list()
	for genotype_label, genotype in modified_genotypes.iterrows():
		if genotype_label in children:
			genotype_frequencies = genotype[genotype > detection_cutoff]
			genotype_children: pandas.Series = modified_genotypes.loc[children[genotype_label]].max()
			genotype_frequencies: pandas.Series = genotype_frequencies - genotype_children
			genotype_frequencies = genotype_frequencies.mask(lambda s: s < detection_cutoff, 0.05)  # 0.05 so there is still a visible slice.
			genotype_frequencies = genotype_frequencies.fillna(0)  # Otherwise plotting the muller diagram will fail.
		else:
			genotype_frequencies = genotype
		genotype_frequencies.name = genotype_label
		children_table.append(genotype_frequencies)
	return pandas.DataFrame(children_table)


def _convert_genotype_table_to_population_table(genotype_table: pandas.DataFrame) -> pandas.DataFrame:
	""" Pivots a genotype table into a long-form table.
		Parameters
		----------
		genotype_table: pandas.DataFrame
			 The genotype table to pivot.

		Returns
		-------
		pandas.DatFrame
			- columns
				`Generation`: int
				`Identity`: str
				`Population`: float
	"""
	table = list()
	for genotype_label, genotype_frequencies in genotype_table.iterrows():
		# Remove timepoints where the value was 0 or less than 0 due to above line.
		for timepoint, frequency in genotype_frequencies.items():
			row = {
				'Identity':   genotype_label,
				'Generation': int(timepoint),
				'Population': frequency * 100
			}
			table.append(row)

	temp_df = pandas.DataFrame(table)
	return temp_df


def _append_genotype_0(population_table: pandas.DataFrame) -> pandas.DataFrame:
	generation_groups = population_table.groupby(by = 'Generation')
	modified_population = population_table
	for generation, group in generation_groups:
		population = group['Population'].sum()
		if population <= 100:
			p = 100 - population
		else:
			p = 0
		modified_population = modified_population.append({'Generation': generation, "Identity": "genotype-0", "Population": p}, ignore_index = True)
	return modified_population


def generate_ggmuller_population_table(mean_genotypes: pandas.DataFrame, edges: pandas.DataFrame, detection_cutoff: float,
		adjust_populations: bool) -> pandas.DataFrame:
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
	adjust_populations:bool
		If true, subtracts any children genotypes from the parent's frequency. Used to make the muller diagrams more intuitive.
	Returns
	-------
	pandas.DataFrame
	"""

	# "Generation", "Identity" and "Population"
	# Adjust populations to account for inheritance.
	# If a child genotype fixed, the parent genotype should be replaced.

	# Use a copy of the dataframe to avoid making changes to the original.
	modified_genotypes: pandas.DataFrame = mean_genotypes.copy(True)
	# In case the genotype table includes genotypes that the edges table does not have.
	modified_genotypes = modified_genotypes[modified_genotypes.index.isin(edges['Identity'])]

	# Generate a list of all genotypes that arise in the background of each genotype.
	# Should ba a dict mapping parent -> list[children]
	children = _compile_parent_linkage(edges)
	if adjust_populations:
		child_df = _subtract_children_from_parent(modified_genotypes, children, detection_cutoff)
	else:
		child_df = modified_genotypes
	temp_df = _convert_genotype_table_to_population_table(child_df)
	population_table = _append_genotype_0(temp_df)

	return population_table


def generate_ggmuller_edges_table(genotype_clusters: Dict[str,List[str]]) -> pandas.DataFrame:
	table = list()
	for identity, background in genotype_clusters.items():
		parent = background[0]
		if parent == identity:
			parent = 'genotype-0'

		row = {
			'Parent':   parent,
			'Identity': identity
		}
		table.append(row)
	table = pandas.DataFrame(table)[['Parent', 'Identity']]  # Reorder columns

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
		pair_df['X'] = calculation.X

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


def generate_trajectory_table(trajectories: pandas.DataFrame, parent_genotypes: Union[pandas.Series, Dict[str, str]],
		info: pandas.DataFrame) -> pandas.DataFrame:
	# Sorts the trajectory table and adds additional columns from the original table.
	trajectories = trajectories[sorted(trajectories.columns, key = lambda s: int(s))]
	if info is not None:
		trajectory_table: pandas.DataFrame = trajectories.copy()
		trajectory_table['genotype'] = [parent_genotypes[k] for k in trajectory_table.index]
		trajectory_table = trajectory_table.join(info).sort_values(by = ['genotype'])
		return trajectory_table


if __name__ == "__main__":
	from pprint import pprint
	from import_data import import_table_from_string

	test_table = import_table_from_string(
		"""
			Parent	Identity
			genotype-0	genotype-1
			genotype-1	genotype-4
			genotype-4	genotype-5
			genotype-0	genotype-2
			genotype-1	genotype-3
		"""
	)
	output = _compile_parent_linkage(test_table)
	pprint(output)
