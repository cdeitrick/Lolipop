from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
# plt.switch_backend('agg')
import pandas
import seaborn


def plot_heatmap(data: pandas.DataFrame, filename: Path):
	font = {
		'size': 20
	}
	use_annotations = len(data) < 20  # So the annotations are actually visible
	matplotlib.rc('font', **font)
	figsize = (20, 20)
	fig, ax = plt.subplots(figsize = figsize)
	seaborn.heatmap(
		data,
		ax = ax,
		annot = use_annotations,
		cmap = 'Reds',
		square = True,
		xticklabels = True,
		yticklabels = True
	)
	ax.tick_params(axis = 'both', which = 'major', labelsize = 24)
	ax.set_ylabel("Trajectory Label", size = 32)
	ax.set_xlabel("Trajectory Label", size = 32)
	ax.set_title("p-values of all mutational trajectories", size = 36)
	plt.tight_layout()
	fig.savefig(str(filename), format = 'svg')
