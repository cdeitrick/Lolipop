import matplotlib.pyplot as plt
import pandas
from muller import varycolor
from dataclasses import dataclass
from typing import List


def get_numeric_columns(columns):
	for i in columns:
		try:
			int(i)
			yield i
		except:
			pass

def get_cluster_color(element, cluster_colors)->str:

	if isinstance(element, (int,float)):
		color = [i for i in cluster_colors if element in i[0]]
	else:
		keys = [int(i) for i in element.split('|')]
		color = [i for i in cluster_colors if set(keys) == set(i[0])]

	return color[0][1]

def plot_genotypes(timeseries: pandas.DataFrame, genotypes: pandas.DataFrame):

	numeric_columns = list(get_numeric_columns(timeseries.columns))
	genotype_colors = [(g, varycolor.generate_random_color()) for g in genotypes]
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
	mean_genotypes = get_mean_genotypes(genotypes, time_series)
	plt.subplot(2, 1, 2)
	plt.ylabel('frequency')
	plt.xlabel('time')
	plt.title('Mean of Genotypes')
	for cluster_id, cluster_timeseries in mean_genotypes.iterrows():
		cluster_color = get_cluster_color(cluster_id, genotype_colors)

		plt.plot(cluster_timeseries.index, cluster_timeseries.values, '-o', color = cluster_color, label = cluster_id)
	plt.legend()
	plt.show()


if __name__ == "__main__":
	import variables
	from time_series_import import import_timeseries
	from get_genotypes import get_genotypes
	from avetrajectories import get_mean_genotypes
	time_series, info = import_timeseries(variables.filename)
	genotypes = get_genotypes(time_series)

	plot_genotypes(time_series, genotypes)
