from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

try:
	from muller.clustering.metrics import DistanceCache
except ModuleNotFoundError:
	from ..clustering.metrics import DistanceCache


def plot_dendrogram(linkage_table: Any, labels, filename: Path):
	# noinspection PyUnusedLocal
	fig, ax = plt.subplots(figsize = (15, 15))
	ax: plt.Axes
	linkage_table = linkage_table[['left', 'right', 'distance', 'observations']].values # Removes extra column
	# plt.figure(figsize = (15, 15))
	ax.set_title('Hierarchical Clustering Dendrogram', size = 40)
	ax.set_xlabel('Trajectory Label', size = 32)
	ax.set_ylabel('Distance', size = 32)
	hierarchy.dendrogram(
		linkage_table,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = labels,
		ax = ax
	)
	ax.tick_params(axis = 'both', labelsize = 20)

	plt.savefig(filename, dpi = 500)
