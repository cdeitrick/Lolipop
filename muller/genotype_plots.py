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
	if genotype_colors is None:
		color_palette = [
			'#bf0000', '#ffa280', '#402b00', '#98b32d', '#3df29d', '#60acbf', '#000f73',
			'#ce3df2', '#e639ac', '#bf6060', '#594943', '#f2ca79', '#4a592d', '#003322',
			'#001b33', '#333366', '#66005f', '#a60058', '#331a1a', '#ff6600', '#d9c7a3',
			'#19bf00', '#468c75', '#266399', '#070033', '#cc99c2', '#ff0044', '#f2beb6',
			'#a65b29', '#8c7000', '#a0cc99', '#00d9ca', '#80c4ff', '#3000b3', '#40303d',
			'#73565e', '#661b00', '#ffaa00', '#ffee00', '#00731f', '#394b4d', '#0066ff',
			'#8f66cc', '#330022', '#59161f'
		]
		genotype_colors = {k:v for k,v in zip(sorted(mutational_genotypes.index), color_palette)}

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
	from pathlib import Path

	try:
		from muller.import_table import import_trajectory_table
		from muller import calculate_genotypes
	except ModuleNotFoundError:
		# noinspection PyUnresolvedReferences
		from import_table import import_trajectory_table
		# noinspection PyUnresolvedReferences
		import calculate_genotypes

	input_filename = Path('/home/cld100/Documents/github/muller_diagrams/Data files/B1_muller_try1.xlsx')

	timepoints, info = import_trajectory_table(input_filename)

	mean_genotypes = calculate_genotypes.workflow(timepoints, calculate_genotypes.GenotypeOptions.from_breakpoints(0.03))
	plot_genotypes(timepoints, mean_genotypes)
