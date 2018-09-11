import matplotlib.pyplot as plt
import pandas

from typing import Union
import random
def generate_random_color()->str:
	r = random.randint(100,255)
	g = random.randint(100,255)
	b = random.randint(100,255)
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

def get_cluster_color(element:Union[int,float, str], cluster_colors)->str:
	"""
		Returns the color of the cluster this genotype belongs to.
	"""
	return generate_random_color()
	if isinstance(element, (int,float)):
		color = [i for i in cluster_colors if element in i[0]]
	else:
		keys = [int(i) for i in element.split('|')]
		color = [i for i in cluster_colors if set(keys) == set(i[0])]

	return color[0][1]

def plot_genotypes(timeseries: pandas.DataFrame, mean_genotypes:pandas.DataFrame):
	"""
		Plots the clustered genotypes.
	Parameters
	----------
	timeseries
	mean_genotypes
	genotypes

	Returns
	-------

	"""
	mean_genotypes.pop('members')
	numeric_columns = list(get_numeric_columns(timeseries.columns))
	genotype_colors = [(g, generate_random_color()) for g in mean_genotypes.index]
	# Plot all trajectories
	plt.subplot(2, 1, 1)

	plt.ylabel('frequency')

	plt.title('All Trajectories')
	for index, trajectory in timeseries.iterrows():

		trajectory_id = trajectory['Trajectory']
		color = get_cluster_color(trajectory_id, genotype_colors)
		trajectory_timeseries = sorted((column, trajectory[column]) for column in numeric_columns)
		x_values, y_values = zip(*trajectory_timeseries)
		plt.plot(x_values, y_values, '-o', color = color, label = trajectory_id)

	#Plot clustered genotypes. Should be same as above, but colored based on genotype cluster.

	# Plot the mean of each cluster
	plt.subplot(2, 1, 2)
	plt.ylabel('frequency')
	plt.xlabel('time')
	plt.title('Mean of Genotypes')
	for cluster_id, cluster_timeseries in mean_genotypes.iterrows():
		cluster_color = get_cluster_color(cluster_id, genotype_colors)

		plt.plot(cluster_timeseries.index, cluster_timeseries.values, '-o', color = cluster_color, label = cluster_id)
	plt.legend()
	plt.show()
	print()


if __name__ == "__main__":
	from pathlib import Path
	try:
		from muller.import_table import import_timeseries
		from muller import get_genotypes
	except ModuleNotFoundError:
		from import_table import import_timeseries
		import get_genotypes

	input_filename = Path('/home/cld100/Documents/github/muller_diagrams/Data files/P1/P1_Muller.xlsx')

	timepoints, info = import_timeseries(input_filename)

	mean_genotypes = get_genotypes.workflow(timepoints)
	plot_genotypes(timepoints, mean_genotypes)
