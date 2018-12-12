import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Dict, List
import random
try:
	from muller.widgets import generate_random_color, get_numeric_columns
except ModuleNotFoundError:
	from widgets import generate_random_color, get_numeric_columns

def plot_genotypes(timeseries: pandas.DataFrame, mutational_genotypes: pandas.DataFrame, filename: Path, genotype_colors:Dict[str,str], parent_genotypes:Dict[str,str]):
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
	Returns
	-------

	"""
	if 'members' in mutational_genotypes:
		genotype_members = mutational_genotypes.pop('members')

	if timeseries is not None:
		numeric_columns = list(get_numeric_columns(timeseries.columns))
	else:
		numeric_columns = list(get_numeric_columns(mutational_genotypes.columns))

	# Plot all trajectories
	grid = plt.GridSpec(2, 3, wspace = 0.4, hspace = 0.3)
	if timeseries is not None:
		plt.subplot(grid[0,0:2])
		plt.ylabel('frequency')
		plt.title('Trajectories')

		for trajectory_id, trajectory in timeseries.iterrows():
			trajectory_genotype_id = parent_genotypes.get(trajectory_id, 'removed')
			color = genotype_colors[trajectory_genotype_id]
			trajectory_timeseries = sorted((column, trajectory[column]) for column in numeric_columns)
			x_values, y_values = zip(*trajectory_timeseries)
			plt.plot(x_values, y_values, '-o', color = color, label = trajectory_id, markersize = 2)

	# Plot clustered muller_genotypes. Should be same as above, but colored based on genotype cluster.
	# Plot the mean of each cluster
	if timeseries is not None:
		plt.subplot(grid[1, :])
	else:
		plt.subplot(grid[:,:])
	plt.ylabel('frequency')
	plt.xlabel('timepoint')
	plt.title('Genotypes')
	for cluster_id, cluster_timeseries in mutational_genotypes.iterrows():
		cluster_color = genotype_colors[cluster_id]
		plt.plot(
			cluster_timeseries.index,
			cluster_timeseries.values,
			'-o', color = cluster_color,
			label = cluster_id,
			markersize = 2
		)
	# noinspection PyUnboundLocalVariable
	if len(cluster_timeseries) > 20:
		bbox = (1, 1.8)
	else:
		bbox = (1, 2)
	legend_font_properties = {'size': 5}
	plt.legend(
		loc = 'right',
		ncol = 2,
		bbox_to_anchor=bbox,
		prop = legend_font_properties,
		title = 'Genotypes',
	)

	if filename:
		plt.savefig(str(filename), dpi = 500, format = 'png')
	else:
		plt.show()



if __name__ == "__main__":
	pass
