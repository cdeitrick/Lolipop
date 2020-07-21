import math
from functools import partial
from typing import *

import pandas
from loguru import logger

Numeric = Union[int, float]

def lag_gens(x: int, values: Iterable[int]) -> int:
	# function to get the generation previous to a specified generation:
	# Implemented in-line so that we can use the population['Genotypes'] within the method.
	candidates = [i for i in values if i < x]
	if candidates:
		return max(candidates)
	else:
		return 0


def dup(iterable: List[Any]) -> List[bool]:
	""" returns a list indicating whether the coresponding index in `iterable` is a duplicated value."""
	return [iterable.count(i) > 1 for i in iterable]


class Vector:
	""" Emulates the way R values can be added to a scalar variable."""

	def __init__(self, value: Optional[Any]):
		self.data: List[Any] = list()
		if value:
			self.data.append(value)

	def __len__(self)->int:
		return len(self.data)

	def __setitem__(self, key: int, value: Any):
		# Inserts `value` at the position `key`
		if key < len(self):
			self.data[key] = value
		elif key - len(self) <= 1:
			self.data.append(value)
		else:
			message = f"Trying to access a value that is not the next position in the array: {key} -> {len(self)}"
			raise ValueError(message)

	def __getitem__(self, item)->Any:
		return self.data[item]


class AdjacencyMatrix:
	def __init__(self, edges: pandas.DataFrame, ages: Optional[pandas.Series] = None):
		self.table_edges = edges

		if ages is not None:
			self.ages = ages
			self.clean_edges = self.generate_clean_edges(ages)
		else:
			self.ages = None
			self.clean_edges = None

	def _validate_point(self, value: int):
		if value not in self.clean_edges['Identity'].values and value not in self.clean_edges['Parent'].values:
			message = f"This is an invalid value: {value}."
			raise ValueError(message)

	def generate_clean_edges(self, lookup: pandas.Series) -> pandas.DataFrame:
		""" This is a more concise implementation of the method used in ggmuller.
			'Replace each genotype name in adjacency matrix with corresponding Age'
		"""
		intermediate_edges = self.table_edges.copy(deep = True)  # To avoid unintented changes.
		intermediate_edges['Identity'] = [lookup.loc[i] for i in self.table_edges['Identity'].values]
		intermediate_edges['Parent'] = [lookup.loc[i] for i in self.table_edges['Parent'].values]
		edges = intermediate_edges.sort_values(by = "Identity")
		self.clean_edges = edges
		self.ages = lookup
		return edges

	def find_start_node(self):
		# Attempts to find the first node in the matrix
		start = self.clean_edges['Parent'].sort_values().iloc[0]  # Just a guess

		while True:
			value = self.move_up(start)
			if value == start: return value
			else:
				start = value

	def move_up(self, identity: int) -> int:
		"""
		Move to parent in adjacency matrix

		Returns the corresponding Parent value.
		When there is no parent (i.e. at the top of the tree), returns the input Identity.
		"""
		# Check if the identity is valid
		self._validate_point(identity)

		parent = self.clean_edges[self.clean_edges["Identity"] == identity].squeeze()  # To make sure it is a series rather than a dataframe
		if parent.empty:
			return identity
		else:
			return parent['Parent']

	def move_down(self, parent: int) -> int:
		self._validate_point(parent)
		daughters: pandas.Series = self.clean_edges[self.clean_edges['Parent'] == parent]['Identity']
		if daughters.empty:
			return parent
		else:
			return daughters.min()

	def move_right(self, identity: int) -> int:
		self._validate_point(identity)

		parent = self.clean_edges[self.clean_edges['Identity'] == identity].squeeze()
		if parent.empty:
			return identity
		parent = parent['Parent']

		siblings = self.clean_edges[self.clean_edges['Parent'] == parent]['Identity'].sort_values().tolist()
		# Find the index of the next value after the one given as the identity.
		try:
			sibling = siblings[siblings.index(identity) + 1]
			return sibling
		except IndexError:
			return identity

	def path_vector(self) -> List[str]:
		n = 0  # Arrays start at zero in python.
		path = self.find_start_node()
		path = Vector(path)  # Since variables in R behave both as scalar and vector values.

		upped = False

		for i in range(10000):  # Prevents endless loops. change to a while loop later
			if not upped:
				for j in range(1000):
					n = n + 1
					path[n] = self.move_down(path[n - 1])
					upped = False
					if path[n] == path[n - 1]: break
			if self.move_right(path[n]) != path[n]:
				n = n + 1
				path[n] = self.move_right(path[n - 1])
				upped = False
			elif self.move_up(path[n]) != path[n]:
				n = n + 1
				path[n] = self.move_up(path[n - 1])
				upped = True

			if path[n] == path[0]: break
			if n > 2 * len(self.clean_edges) + 2:
				message = "Error: stuck in a loop"
				raise ValueError(message)
		# TODO add the final check. Not included yet to save time.from
		if len(path) != (2 * len(self.clean_edges) + 2):
			message = "Error: the adjacency matrix seems to be bipartite"
			raise ValueError(message)

		# Map the path of ages to the corresponding genotype
		path = [self.age_to_genotype(i) for i in path]

		return path

	def age_to_genotype(self, age: int) -> str:
		""" Returns the genotype represented by the given age."""
		if self.ages is None:
			message = f"The lookup table has not been provided yet."
			raise ValueError(message)
		counts = self.ages.value_counts()
		if any(i > 1 for i in counts.values):
			message = "Two or more genotypes were assigned the same age."
			raise ValueError(message)

		reverse_lookup = {v: k for k, v in self.ages.items()}

		return reverse_lookup[age]


class GenerateMullerDataFrame:
	""" implements the get_Muller_df function from ggmuller. The original fails when there are a large number of sampled timepoints."""

	def __init__(self):
		self.time_column = 'Generation'
		self.identity_column = 'Identity'
		self.population_column = 'Population'

		# Initial population size. Will probably never be changed.
		self.initial_population_size = 0
		self.initial_start_positions = 0.5
		self.detection_limit = 0  # Would normally use 0.05 * 100, but the original script uses 0.



	# Add the sub-workflows that will be used.
	@staticmethod
	def _adjust_population_table(population: pandas.DataFrame, first_generation: pandas.DataFrame, start_positions: float) -> pandas.DataFrame:
		# copy all rows for generations at which new genotypes appear
		new_rows = population[population['Generation'].isin(first_generation['start_time'].values)]
		prev_rows = population[population['Generation'].isin(first_generation['previous_time'].values)]
		prev_rows.index = [i + 1 for i in prev_rows.index]

		# adjust generations of copied rows
		partial_lag_gens = partial(lag_gens, values = population['Generation'].unique())
		adjusted_generation:List[float] = new_rows['Generation'] - start_positions * (new_rows['Generation'] - new_rows['Generation'].apply(partial_lag_gens))

		adjusted_population = ((1 - start_positions) * new_rows['Population']) + (start_positions * prev_rows['Population'])

		new_rows['Generation'] = adjusted_generation
		new_rows['Population'] = adjusted_population

		return new_rows.reset_index(drop = True) # To make the tale more consistent

	@staticmethod
	def _clean_population_table(population: pandas.DataFrame) -> pandas.DataFrame:
		""" Implements the cleaning steps from the ggmuller package. Since these are already done with the other scripts, this returns
			The input table.
		"""
		return population


	def _find_start_positions(self, series: pandas.Series, start_positions: float) -> float:
		all_generations = series.unique()
		difference_time_diff = series.diff().abs().min()  # .diff() is reversed from r, so use abs to correct it.
		difference_time_absolute = max(all_generations) - min(all_generations)
		delta = min(1E-2 * difference_time_diff, 1E-4 * difference_time_absolute)

		start_positions = max(start_positions, delta)
		start_positions = min(start_positions, 1 - delta)

		# I think the ggmuller scripts assumed there would be relatively small time intervals, such as `35 days` or `120` minutes
		# The LTEE has several thousand, so the above code results in a negative value. For now replace negative values with the original
		# values of `start_positions`.
		start_positions = max(start_positions, self.initial_start_positions)

		return start_positions

	def _get_initial_generations(self, population: pandas.DataFrame) -> pandas.DataFrame:
		""" Maps each genotype to both the first timepoint it appears as well as the previous timepoint."""
		genotype_groups = population.groupby(by = "Identity")
		table = list()
		for identity, group in genotype_groups:
			detected = group[group[self.population_column] > self.detection_limit]
			detected_timepoint = detected[self.time_column].min()

			# Also need to find the previous timepoint.

			previous_timepoints = [i for i in group[self.time_column].values if i < detected_timepoint]
			if previous_timepoints:  # Check if there are any previous timepoints.
				previous_timepoint = max(previous_timepoints)
			else:
				previous_timepoint = detected_timepoint  # Not sure if this is the way the original script handled this edge case.

			row = {
				self.identity_column: identity,
				'start_time':         detected_timepoint,
				'previous_time':      previous_timepoint
			}
			table.append(row)

		df = pandas.DataFrame(table)

		df = df[[self.identity_column, 'start_time', 'previous_time']]  # sort columns to match ggmuller

		# Remove `generation-0` to match the ggmuller script
		df = df[df['start_time'] != 0]
		return df.reset_index(drop = True) # The index isn't needed, so reset it to make it more consistent with the other tables.

	@staticmethod
	def expand(left_values: Iterable[Any], right_values: Iterable[Any]) -> pandas.DataFrame:
		""" Generates a two-column table which pairs every left value with every right value.
			The resulting table will be len(left) * len(right) rows.
		"""
		added_rows = list()
		for left in left_values:
			for right in right_values:
				added_rows.append({'Left': left, 'Right': right})
		added_rows = pandas.DataFrame(added_rows)

		return added_rows

	@staticmethod
	def add_semi_frequencies(population: pandas.DataFrame) -> pandas.DataFrame:
		logger.warning("patch")
		print(population.shape)
		print(population.to_string())

		population = population.dropna()
		frequencies = population.groupby(by = "Generation").apply(lambda s: .5 * (s['Population'] / s['Population'].sum()))
		print()
		print(frequencies.shape)
		print(frequencies.to_string())
		population['Frequency'] = frequencies.fillna(0).values
		population['Population'] = population['Population'] / 2  # Because of the duplication

		return population

	def apply_add_start_points(self, population: pandas.DataFrame, start_positions: float = 0.5) -> pandas.DataFrame:
		start_positions = self._find_start_positions(population['Generation'], start_positions)
		# set small initial population size
		first_generation = self._get_initial_generations(population)

		# if all genotypes appear at the first time point then don't make any changes:
		if len(first_generation['start_time'].unique()) == 1:
			return population

		# Adjust the generations and populations of population
		new_rows = self._adjust_population_table(population, first_generation, start_positions)

		# Add the new columns to the original population dataframe
		pop_df = pandas.concat([population, new_rows])
		pop_df = pop_df.sort_values(by = ["Generation", "Identity"]).reset_index(drop = True)

		# adjust generations in reference list
		first_generation['Generation'] = first_generation["start_time"] - start_positions * (
				first_generation['start_time'] - first_generation['previous_time'])
		# Set small initial populations in reference list
		first_generation['Population2'] = self.initial_population_size

		# replace initial populations in the dataframe with values from reference list
		pop_df = pop_df.merge(first_generation, how = 'outer')
		# Replace the population of rows where `Population2` is defined.
		pop_df['Population'] = [(i if math.isnan(j) else j) for i, j in zip(pop_df['Population'].values, pop_df['Population2'].values)]

		return pop_df[['Generation', 'Identity', 'Population']]

	def get_genotype_ages(self, table: pandas.DataFrame) -> pandas.Series:
		# Apply the filtering step before the groupby step.

		df = table[(table['Population'] > self.detection_limit) | (table['Generation'] == table['Generation'].max())]

		groups = df.groupby(by = 'Identity')

		first = groups.first()
		df = first.sort_values(by = ["Generation", "Identity"], ascending = [True, True])
		df = df.reset_index()
		df['Age'] = list(range(1, len(df) + 1))
		series = df.set_index("Identity")['Age']
		return series

	def reorder_by_vector(self, muller_df: pandas.DataFrame, path_vector: List[str]):
		""" Orders the muller dataframe so that each series is plotted in the correct order."""

		# Add a unique id column to the vector.
		# Identify duplicates
		# Using pandsa.Series.value_counts() might be easier, but would require makeing a Series object.

		vector_df = pandas.DataFrame({
			"vector":           path_vector,
			'uniqueGenotypeId': self._reorder_by_vector_generate_unique_ids_genotype(path_vector)
		})

		# Now start modifying `muller_df`
		# Take the vector with the unique genotype id's and multiply it so that it is the same length as `muller_df`
		number_of_unique_timepoints = len(muller_df['Generation'].unique())
		vector: List[str] = vector_df['uniqueGenotypeId'].tolist() * number_of_unique_timepoints
		# Assume that the `muller_df` is already sorted by generation.
		# Basically, each repeated sequence in `vector` corresponds to a single generation.
		unique_id_generation = ["{}_{}".format(i, g if g != int(g) else int(g)) for i, g in zip(vector, muller_df['Generation'].tolist())]

		# Index the dataframe by the unique generation id so that .loc[] can be used to reorder it.

		# Name it `Unique_id` to match the ggmuller scripts.
		muller_df['Unique_id'] = self._reorder_by_vector_generate_unique_ids_generation(muller_df)
		# To match the scripts
		muller_df['Group_id'] = muller_df['Unique_id'].apply(lambda s: s.split('_')[0])
		muller_df_indexed = muller_df.set_index('Unique_id')
		muller_df_sorted = muller_df_indexed.loc[unique_id_generation]

		return muller_df_sorted.reset_index()

	@staticmethod
	def _reorder_by_vector_generate_unique_ids_generation(df: pandas.DataFrame) -> List[str]:
		seen = set()
		uids = list()
		for index, row in df.iterrows():
			name = row['Identity']
			generation = row['Generation']
			# Check if the `generation` value can be safely converted to int to match ggmuller.
			if generation - int(generation) == 0:
				generation = int(generation)
			key = (name, generation)
			if key in seen:
				name += "a"
			else:
				seen.add(key)
			uid = f"{name}_{generation}"
			uids.append(uid)
		return uids

	@staticmethod
	def _reorder_by_vector_generate_unique_ids_genotype(iterable: List[str]) -> List[str]:
		"""
		Since each series has been split in two, there should be exactly one duplicate for each value.
		Since each value needs to have a unique id, append a 'a' to the second value.
		Here is a quick and dirty implementation:
		"""
		seen = set()
		unique_ids = list()
		for element in iterable:
			if element in seen:
				unique_ids.append(element + 'a')
			else:
				unique_ids.append(element)
				seen.add(element)
		return unique_ids

	def correct_population_values(self, population: pandas.DataFrame) -> pandas.DataFrame:
		""" Adjusts the values in the population dataframe to make sure it can be plotted correctly."""
		# add rows to population to ensure genotype starting positions are plotted correctly.
		population = self.apply_add_start_points(population, self.initial_start_positions)
		# Construct a dataframe with the "Age" of each genotype.
		population_sorted = population.sort_values(by = ["Generation", "Population"], ascending = [True, False])
		# Add semi-frequencies
		frequencies = self.add_semi_frequencies(population_sorted)

		# duplicate the rows on the dataframe since the muller plots need to plot the top and bottom halves of each series individually.
		muller_df = pandas.concat([frequencies, frequencies]).sort_values(by = 'Generation')

		return muller_df

	def run(self, edges: pandas.Series, population: pandas.DataFrame) -> pandas.DataFrame:
		# Test if the length of the population table makes sense given the population and number of timepoints.
		edges = edges.reset_index()
		population = self._clean_population_table(population)
		# Construct a dataframe with the "Age" of each genotype.
		population_sorted = population.sort_values(by = ["Generation", "Population"], ascending = [True, False])

		# Generate a dataframe with the adjusted generations and populations.
		muller_df = self.correct_population_values(population)

		ages = self.get_genotype_ages(population_sorted)
		clean_edges = AdjacencyMatrix(edges, ages)
		path_vector: List[str] = clean_edges.path_vector()
		path_vector: List[str] = path_vector[::-1]  # Convert to list and reverse it to match ggmuller scripts.
		muller_ordered = self.reorder_by_vector(muller_df, path_vector)
		# Replace each age in the path vector with the corresponding genotype name
		# Rearrange the columns so pandas.testing.assert_frame_equal doesn't complain.
		muller_ordered = muller_ordered[['Generation', 'Identity', 'Population', 'Frequency', 'Group_id', 'Unique_id']]
		return muller_ordered
