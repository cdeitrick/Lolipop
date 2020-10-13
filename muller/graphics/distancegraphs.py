from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import seaborn


# plt.style.use("/home/cld100/Documents/sandbox/matplotlibrc")


def generate_distance_plot(distances: List[float], similarity_cutoff: float, filename: Optional[Path] = None, ax: plt.Axes = None):
	""" Shows the spread of computed pairwise distances as a histogram."""
	if not ax:
		fig, ax = plt.subplots(figsize = (12, 10))

	seaborn.distplot(distances, ax = ax, kde = False, rug = True, bins = 20)
	ax.axvline(similarity_cutoff, color = 'red')

	ax.set_title("Pairwise distances between each pair of trajectories")
	ax.set_xlabel("Distance")
	ax.set_ylabel("Count")
	ax.set_xlim(0, max(distances))
	plt.tight_layout()
	if filename:
		plt.savefig(filename)
	else:
		plt.show()
