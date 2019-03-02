from pathlib import Path
from typing import Any, Dict

import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

try:
	from clustering.metrics import PairwiseCalculationCache
except ModuleNotFoundError:
	from ..clustering.metrics import PairwiseCalculationCache


def get_dendrogram_colors(Z, trajectory_colors: Dict[str, str]) -> Dict[str, str]:
	# see question for code prior to "color mapping"

	# Color mapping
	trajectory_colors["0"] = "#FF0000"
	dflt_col = "#808080"  # Unclustered gray
	# notes:
	# * rows in Z correspond to "inverted U" links that connect clusters
	# * rows are ordered by increasing distance
	# * if the colors of the connected clusters match, use that color for link
	index_to_trajectory = lambda s: str(s)

	link_cols: Dict[str, str] = {}
	for cluster_index, (left_label, right_label) in enumerate(Z[:, :2].astype(int)):
		get_color = lambda s: link_cols[s] if s > len(Z) else trajectory_colors[index_to_trajectory(s)]
		color_1 = get_color(left_label)
		color_2 = get_color(right_label)
		link_cols[cluster_index + 1 + len(Z)] = color_1 if color_1 == color_2 else dflt_col

	return link_cols


def plot_dendrogram(Z: Any, pair_array: PairwiseCalculationCache, filename: Path):
	fig, ax = plt.subplots(figsize = (15, 15))
	#plt.figure(figsize = (15, 15))
	ax.set_title('Hierarchical Clustering Dendrogram', size = 30)
	ax.set_xlabel('Trajectory Label', size = 20)
	ax.set_ylabel('Distance', size = 20)
	hierarchy.dendrogram(
		Z,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = pair_array.squareform().index,
		ax = ax
	)

	plt.savefig(filename, dpi = 500)
