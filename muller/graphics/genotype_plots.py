import matplotlib.pyplot as plt
from matplotlib.figure import Axes  # For autocomplete

plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Dict, Optional

try:
	from muller.widgets import get_numeric_columns
except ModuleNotFoundError:
	from ..widgets import get_numeric_columns
BASE_COLOR = "#333333"


def plot_timeseries(timeseries: pandas.DataFrame, palette: Dict[str, str], ax: Optional[Axes] = None, filename: Optional[Path] = None):
	if ax is None:
		fig, ax = plt.subplots(figsize = (10, 5))
	ax.set_ylabel('frequency')
	plot_title = 'Genotypes' if 'genotype' in timeseries.index[0] else 'Trajectories'
	ax.set_title(plot_title)
	ax.set_ylim(0, 1.01)

	numeric_columns = list(get_numeric_columns(timeseries.columns))

	for trajectory_id, trajectory in timeseries.iterrows():
		color = palette.get(trajectory_id, BASE_COLOR)
		trajectory_timeseries = sorted((column, trajectory[column]) for column in numeric_columns)
		x_values, y_values = zip(*trajectory_timeseries)
		ax.plot(x_values, y_values, '-o', color = color, label = trajectory_id, markersize = 2)
	if filename:
		plt.tight_layout()
		plt.savefig(str(filename), dpi = 200)
	return ax


# noinspection PyTypeChecker
def plot_genotypes(timeseries: pandas.DataFrame, mutational_genotypes: pandas.DataFrame, filename: Path, genotype_colors: Dict[str, str],
		trajectory_palette: Dict[str, str]):
	"""
		Plots the clustered muller_genotypes.
	Parameters
	----------
	timeseries: pandas.DataFrame
		The table with each trajectory used to generate the muller_genotypes.
	filename: Optional[Path]
		Path to save the genotype plot.
	mutational_genotypes: pandas.DataFrame
	genotype_colors: Dict[str,str]
		A mapping of muller_genotypes to their corresponding colors.
	trajectory_palette: Dict[str,str]
		Maps trajectories to their respective genotype color.
	Returns
	-------

	"""
	if 'members' in mutational_genotypes:
		mutational_genotypes.pop('members')

	# Plot all trajectories
	grid = plt.GridSpec(2, 3, wspace = 0.4, hspace = 0.3)
	if timeseries is not None:
		taxes: Axes = plt.subplot(grid[0, 0:2])
		plot_timeseries(timeseries, trajectory_palette, taxes)

	# Plot clustered muller_genotypes. Should be same as above, but colored based on genotype cluster.
	# Plot the mean of each cluster
	if timeseries is not None:
		gaxes: Axes = plt.subplot(grid[1, :])
	else:
		gaxes: Axes = plt.subplot(grid[:, :])

	plot_timeseries(mutational_genotypes, genotype_colors, gaxes)
	# noinspection PyUnboundLocalVariable
	if len(mutational_genotypes) > 20:
		bbox = (1, 1.8)
	else:
		bbox = (1, 2)
	legend_font_properties = {'size': 5}
	plt.legend(
		loc = 'right',
		ncol = 2,
		bbox_to_anchor = bbox,
		prop = legend_font_properties,
		title = 'Genotypes',
	)

	if filename:
		plt.savefig(str(filename), dpi = 500, format = 'png')
	else:
		plt.show()
