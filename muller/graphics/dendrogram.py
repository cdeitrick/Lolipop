from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

try:
	from clustering.metrics import PairwiseCalculationCache
except ModuleNotFoundError:
	from ..clustering.metrics import PairwiseCalculationCache


def plot_dendrogram(linkage_table: Any, pair_array: PairwiseCalculationCache, filename: Path):
	fig, ax = plt.subplots(figsize = (15, 15))
	# plt.figure(figsize = (15, 15))
	ax.set_title('Hierarchical Clustering Dendrogram', size = 30)
	ax.set_xlabel('Trajectory Label', size = 20)
	ax.set_ylabel('Distance', size = 20)
	hierarchy.dendrogram(
		linkage_table,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = pair_array.squareform().index,
		ax = ax
	)

	plt.savefig(filename, dpi = 500)
