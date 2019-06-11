from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

try:
	from muller.clustering.metrics import PairwiseCalculationCache
except ModuleNotFoundError:
	from ..clustering.metrics import PairwiseCalculationCache


def plot_dendrogram(linkage_table: Any, pair_array: PairwiseCalculationCache, filename: Path):
	# noinspection PyUnusedLocal
	fig, ax = plt.subplots(figsize = (15, 15))
	ax: plt.Axes
	# plt.figure(figsize = (15, 15))
	ax.set_title('Hierarchical Clustering Dendrogram', size = 40)
	ax.set_xlabel('Trajectory Label', size = 32)
	ax.set_ylabel('Distance', size = 32)
	hierarchy.dendrogram(
		linkage_table,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = pair_array.squareform().index,
		ax = ax
	)
	for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
		label.set_fontsize(20)

	plt.savefig(filename, dpi = 500)
