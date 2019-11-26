from typing import Dict, List, Union

import pandas


class GGMuller:
	""" Consolidates the functions used to save the ggmuller tables.
		Parameters
		----------
		cutoff_detection: float
		adjust_populations:bool
			If true, subtracts any children genotypes from the parent's frequency. Used to make the muller diagrams more intuitive.
			This will make the muller plot show each genotype that reaches 100% as actually taking up the full vertical space in the graphic,
			Rather than as 50%/50% alongside any child genotypes.
	"""

	def __init__(self, cutoff_detection: float, adjust_populations: bool = True):
		self.cutoff_detection = cutoff_detection
		self.adjust_populations = adjust_populations

		# This is the amount to preserve when adjusting the populations frequencies to account for child genotypes.
		# This will allow the parent genotypes to essentially "outline" the child genotypes.
		self.visible_slice = 0.01
		self.ancestral_genotype_label = 'genotype-0'  # The name to assign to the ancestral genotype.

	@staticmethod
	def _compile_parent_linkage(edges: pandas.Series) -> Dict[str, List[str]]:
		""" Maps a genotype to a list of all genotypes that inherit from it"""
		children = dict()
		for identity, parent in edges.items():
			if parent not in children:
				children[parent] = list()
			children[parent].append(identity)
		return children

	@staticmethod
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

	def _subtract_children_from_parent(self, modified_genotypes: pandas.DataFrame, children: Dict[str, List[str]]) -> pandas.DataFrame:
		""" Reduces the observed frequency of parent genotypes to allow child genotypes to be visible at the correct vertical abundance."""

		children_table = list()
		for genotype_label, genotype in modified_genotypes.iterrows():
			if genotype_label in children:
				genotype_frequencies = genotype[genotype > self.cutoff_detection]
				genotype_children: pandas.Series = modified_genotypes.loc[children[genotype_label]].max()
				genotype_frequencies: pandas.Series = genotype_frequencies - genotype_children
				genotype_frequencies = genotype_frequencies.mask(lambda s: s < self.cutoff_detection, self.visible_slice)
				genotype_frequencies = genotype_frequencies.fillna(0)  # Otherwise plotting the muller diagram will fail.
			else:
				genotype_frequencies = genotype
			genotype_frequencies.name = genotype_label
			children_table.append(genotype_frequencies)
		return pandas.DataFrame(children_table)

	def add_ancestral_genotype(self, population_table: pandas.DataFrame) -> pandas.DataFrame:
		""" Adds the ancestral genotype to the graphic. This is based on the observed frequencies at each timepoint and whether they sum to 100%."""
		generation_groups = population_table.groupby(by = 'Generation')
		modified_population = population_table
		for generation, group in generation_groups:
			# Check if the observed genotypes collectively sum to 100%.
			population = group['Population'].sum()
			if population > 100:
				p = 0
			else:
				p = 100 - population
			current_timepoint = {
				'Generation': generation,
				"Identity":   self.ancestral_genotype_label,
				"Population": p
			}
			# Since `modified_population` is a dataframe, need to manually tell it to ignore the index.
			modified_population = modified_population.append(current_timepoint, ignore_index = True)
		return modified_population

	def generate_ggmuller_population_table(self, edges: pandas.Series, mean_genotypes: pandas.DataFrame) -> pandas.DataFrame:
		"""
			Converts the genotype frequencies to a population table suitable for ggmuller.
		Parameters
		----------
		mean_genotypes: pandas.DataFrame
			Contains the mean frequency for each genotype at each timepoint.
		edges: pandas.Series
			A Mapping from child genotype to parent genotype.

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
		modified_genotypes = modified_genotypes[modified_genotypes.index.isin(edges.index)]

		# Generate a list of all genotypes that arise in the background of each genotype.
		# Should ba a dict mapping parent -> list[children]
		children = self._compile_parent_linkage(edges)
		if self.adjust_populations:
			child_df = self._subtract_children_from_parent(modified_genotypes, children)
		else:
			child_df = modified_genotypes
		temp_df = self._convert_genotype_table_to_population_table(child_df)
		population_table = self.add_ancestral_genotype(temp_df)

		return population_table


def generate_trajectory_table(trajectories: pandas.DataFrame, parent_genotypes: Union[pandas.Series, Dict[str, str]],
		info: pandas.DataFrame) -> pandas.DataFrame:
	# Sorts the trajectory table and adds additional columns from the original table.
	trajectories = trajectories[sorted(trajectories.columns, key = lambda s: int(s))]
	if info is not None:
		trajectory_table: pandas.DataFrame = trajectories.copy()
		trajectory_table['genotype'] = [parent_genotypes.get(k) for k in trajectory_table.index]
		trajectory_table = trajectory_table.join(info).sort_values(by = ['genotype'])
		return trajectory_table
