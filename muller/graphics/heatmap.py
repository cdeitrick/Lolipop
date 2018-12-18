
import matplotlib.pyplot as plt
import matplotlib
#plt.switch_backend('agg')
import pandas
from pathlib import Path
try:
	import seaborn
except ModuleNotFoundError:
	# Seaborn is not installed. Skip plotting
	seaborn = None
def plot_heatmap(data:pandas.DataFrame, filename:Path):
	font = {
		'size':   20
	}

	matplotlib.rc('font', **font)
	figsize = (20, 20)
	fig, ax = plt.subplots(figsize = figsize)
	ax.set_ylabel("Trajectory Label")
	ax.set_xlabel("Trajectory Label")
	ax.set_title("p-values of all mutational trajectories")
	if seaborn is not None:
		seaborn.heatmap(data, ax = ax, annot = True)
		fig.savefig(str(filename), format = 'png')