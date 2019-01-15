import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from typing import Any
from pathlib import Path

def plot_dendrogram(Z:Any, pair_array:Any, filename: Path):
	plt.figure(figsize = (15, 15))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')
	hierarchy.dendrogram(
		Z,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = pair_array.squareform('X').index
	)

	plt.savefig(filename)