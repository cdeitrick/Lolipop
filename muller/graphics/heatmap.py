import seaborn
import matplotlib.pyplot as plt
import pandas
from pathlib import Path
def plot_heatmap(data:pandas.DataFrame, filename:Path):
	figsize = (20, 20)
	fig, ax = plt.subplots(figsize = figsize)
	ax.set_ylabel("Trajectory Label")
	ax.set_xlabel("Trajectory Label")
	ax.set_title("p-values of all mutational trajectories")
	seaborn.heatmap(data, ax = ax)
	fig.savefig(str(filename), format = 'png')