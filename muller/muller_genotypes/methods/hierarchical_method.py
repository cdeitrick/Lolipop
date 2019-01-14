import itertools
from typing import List

import pandas
from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller_genotypes.metrics.similarity import PairwiseArrayType
except ModuleNotFoundError:
	from muller_genotypes.metrics.similarity import PairwiseArrayType


def convert_dictionary_to_squareform(pairs: PairwiseArrayType) -> pandas.DataFrame:
	""" Converts a dictionary with all pairwise values for a set of points into a square matrix representation."""
	keys = sorted(set(itertools.chain.from_iterable(pairs.keys())))
	_square_map = dict()
	for left in keys:
		series = dict()
		for right in keys:
			value = pairs.get((left, right), 0)
			series[right] = value

		_square_map[left] = series

	return pandas.DataFrame(_square_map)


def hierarchical_method(pair_array: PairwiseArrayType, similarity_cutoff: float) -> List[List[str]]:
	numerical_pairs = {k: v.X for k, v in pair_array.items()}
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
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'monocrit', monocrit = MR)

	cluster_map = {i: list() for i in clusters}
	for i, j in zip(clusters, squaremap.index):
		cluster_map[i].append(j)
	return list(cluster_map.values())
