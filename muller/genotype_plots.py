import matplotlib.pyplot as plt
import pandas
from pathlib import Path
from typing import Dict, List
import random

#plt.style.use('classic')

def generate_random_color() -> str:
	r = random.randint(100, 255)
	g = random.randint(100, 255)
	b = random.randint(100, 255)
	color = "#{:>02X}{:>02X}{:>02X}".format(r, g, b)
	return color


def get_numeric_columns(columns):
	"""
		Return a list of all columns which can be converted to a number. These columns represent the timepoints for
		each trajectory. """
	for i in columns:
		try:
			int(i)
			yield i
		except ValueError:
			pass


def find_parent_genotype(trajectory_id: int, genotypes: Dict[str, List[int]]):
	for k, v in genotypes.items():
		if trajectory_id in v:
			return k
	return None


def plot_genotypes(timeseries: pandas.DataFrame, mutational_genotypes: pandas.DataFrame, filename: Path = None, genotype_colors:Dict[str,str]=None):
	"""
		Plots the clustered genotypes.
	Parameters
	----------
	timeseries: pandas.DataFrame
		The table with each trajectory used to generate the genotypes.
	filename: Optional[Path]
		Path to save the genotype plot.
	mutational_genotypes: pandas.DataFrame
	genotype_colors: Dict[str,str]
		A mapping of genotypes to their corresponding colors.
	Returns
	-------

	"""
	genotype_members = mutational_genotypes.pop('members')
	genotype_members = {k: [j for j in v.split('|')] for k, v in sorted(genotype_members.items())}
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
			trajectory_genotype_id = find_parent_genotype(trajectory_id, genotype_members)
			color = genotype_colors[trajectory_genotype_id]
			trajectory_timeseries = sorted((column, trajectory[column]) for column in numeric_columns)
			x_values, y_values = zip(*trajectory_timeseries)
			plt.plot(x_values, y_values, '-o', color = color, label = trajectory_id, markersize = 2)

	# Plot clustered genotypes. Should be same as above, but colored based on genotype cluster.

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
