import itertools
from typing import List, Optional

import pandas
from scipy.cluster import hierarchy
from scipy.spatial import distance
try:
	from muller.muller_genotypes.difference import unlink_unrelated_trajectories
	from muller.muller_genotypes.similarity import PairwiseArrayType
except ModuleNotFoundError:
	from .difference import unlink_unrelated_trajectories
	from .similarity import PairwiseArrayType


def _find_genotype_from_trajectory(element: str, all_genotypes: List[List[str]]) -> Optional[List[str]]:
	"""
		Finds the genotype that contains the trajectory.
	Parameters
	----------
	element: str
		The trajectory id.
	all_genotypes: List[List[int]]
		All muller_genotypes that have been calculated.
	Returns
	-------
		The genotype (in the form of a list of trajectory ids) containing the given trajectory id.
		If the trajectory is not contained in any muller_genotypes, returns None
	"""
	candidates = [i for i in all_genotypes if element in i]
	try:
		value = candidates[0]
	except IndexError:
		value = None

	return value


def _group_trajectories_into_genotypes(pairs: PairwiseArrayType, relative_cutoff: float) -> List[List[str]]:
	"""
		Clusters all trajectories into related muller_genotypes.
		By default the first trajectory makes a genotype category
	Parameters
	----------
	pairs: PairwiseArrayType
		A dictionary mapping pairs to p-values.
	relative_cutoff: float; default 0.05
		The cutoff indicating two trajectories are related.

	Returns
	-------
		A list of all muller_genotypes (each of which is a list of ints)
	"""
	if pairs:
		genotype_candidates = [[min(pairs.keys())[0]]]  # by default the first trajectory forms the first genotype.
	else:
		genotype_candidates = []
	seen = set()
	for (left, right), p_value in pairs.items():
		# ignore pairs that have already been sorted into a genotype.
		if (left, right) in seen or (right, left) in seen:
			continue
		seen.add((left, right))
		# are the genotypes related?
		if p_value > relative_cutoff:
			# Check if any of the trajectories are already listed in genotypes.
			# These will return None if no genotype is found.
			genotype_left = _find_genotype_from_trajectory(left, genotype_candidates)
			genotype_right = _find_genotype_from_trajectory(right, genotype_candidates)

			if genotype_left and genotype_right:
				# they are listed under two different genotypes. Combine them.
				if genotype_left != genotype_right:
					genotype_left += genotype_right
					# Remove the redundant genotype.
					genotype_candidates.remove(genotype_right)

			elif genotype_left:
				genotype_left.append(right)
			elif genotype_right:
				genotype_right.append(left)
			else:
				# Neither element is listed. Create a new genotype
				genotype_candidates.append([left, right])
	return genotype_candidates


def matlab_method(timeseries: pandas.DataFrame, pair_array: PairwiseArrayType, similarity_breakpoint: float, difference_breakpoint: float) -> List[
	List[str]]:
	"""
		Clusters trajectories into muller_genotypes.
	Parameters
	----------
	timeseries: pandas.DataFrame
		The timeseries output from import_timeseries.
			- Index -> str
				Names unique to each trajectory.
			- Columns -> int
				The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	pair_array: PairwiseArrayType
	similarity_breakpoint: float
	difference_breakpoint: float
	Returns
	-------
	List[List[str]]
		A list of the muller_genotypes, where each genotype is itself a list of the name of each member trajectory.
		Ex. [
			[A1, A2, A3],
			[B1, B2],
			[C1]
		]
	"""

	# Trajectories represent the population frequencies at each timepoint
	# Each row represents a single timepoint, each column represents a mutation.
	numerical_array = {k:v.pvalue for k,v in pair_array.items()}
	population_genotypes = _group_trajectories_into_genotypes(numerical_array, similarity_breakpoint)

	# at the end, look at all trajectories that are not listed and
	# append them as their own category.
	# List of all trajectories
	flattened_genotypes = list(itertools.chain.from_iterable(population_genotypes))
	# Any missing trajectories. Will probably be empty
	other_trajectories = [i for i in timeseries.index if i not in flattened_genotypes]
	population_genotypes.append(other_trajectories)

	# finally, for each genotype, make sure each trajectory pair has some
	# non-trivial linkage (say, >0.0005). this avoids falsely linking together
	# trajectories. if not, divide the offending trajectories into two
	# camps, and sort the rest according to which one they are more
	# closely linked to. repeat until everything is linked.
	# pprint(pair_array)

	while True:
		starting_size_of_the_genotype_array = len(population_genotypes)
		population_genotypes = unlink_unrelated_trajectories(population_genotypes[:], numerical_array, difference_breakpoint)
		if len(population_genotypes) == starting_size_of_the_genotype_array:
			break

	# all_genotypes[population_id] = population_genotypes
	# TODO Fix so that it returns a genotype for each population
	return [i for i in population_genotypes if i]  # Only return non-empty lists.


def convert_dictionary_to_squareform(pairs: PairwiseArrayType) -> pandas.DataFrame:
	keys = sorted(set(itertools.chain.from_iterable(pairs.keys())))
	_square_map = dict()
	for left in keys:
		series = dict()
		for right in keys:
			value = pairs.get((left, right), 0)
			series[right] = value

		_square_map[left] = series

	return pandas.DataFrame(_square_map)


def hierarchical_method(pair_array: PairwiseArrayType, similarity_cutoff:float) -> List[List[str]]:

	numerical_pairs = {k:v.X for k, v in pair_array.items()}
	squaremap = convert_dictionary_to_squareform(numerical_pairs)
	condensed_squaremap = distance.squareform(squaremap.values)

	Z = hierarchy.linkage(condensed_squaremap)

	import matplotlib.pyplot as plt
	plt.figure(figsize = (10, 10))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')

	hierarchy.dendrogram(
		Z,
		color_threshold = similarity_cutoff,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = squaremap.index
	)
	plt.savefig("clusterimage.png")
	if False:
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'distance')
	else:
		inconsistent = hierarchy.inconsistent(Z)
		MR = hierarchy.maxRstat(Z, inconsistent, 1)
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion='monocrit', monocrit=MR)
	cluster_map = {i: list() for i in clusters}
	for i, j in zip(clusters, squaremap.index):
		cluster_map[i].append(j)
	return list(cluster_map.values())


if __name__ == "__main__":
	from pprint import pprint
	string = """
		Trajectory	0	1	2	3	4	5
		trajectory-A1	0	0	0	0.1	0.5	0.5
		trajectory-B1	0	0.1	0.15	0.03	0	0
		trajectory-A2	0	0	0	0.06	0.35	0.4
		trajectory-C1	0	0	0	0.3	0.7	1
		trajectory-A3	0	0	0	0	0.45	0.5
		trajectory-B2	0	0.07	0.1	0.02	0.01	0
	"""
	from import_table import import_table_from_string
	data = import_table_from_string(string, index = 'Trajectory')
	from muller.muller_genotypes import similarity
	pair_array = similarity.calculate_pairwise_trajectory_similarity(data, 0.03, 0.97)
	result = hierarchical_method(pair_array)
	pprint(result)

