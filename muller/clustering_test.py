# needed imports

# needed imports
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from muller.muller_genotypes.similarity import calculate_p_value
from muller.import_table import import_table_from_string
from functools import partial
if __name__ == "__main__":
	string = """
		Trajectory	0	1	2	3	4	5
		trajectory-A1	0	0	0	0.1	0.5	0.5
		trajectory-B1	0	0.1	0.15	0.03	0	0
		trajectory-A2	0	0	0	0.06	0.35	0.4
		trajectory-C1	0	0	0	0.3	0.7	1
		trajectory-A3	0	0	0	0	0.45	0.5
		trajectory-B2	0	0.07	0.1	0.02	0.01	0
	"""
	matrix_string = """
	Index	trajectory-C1	trajectory-B1	trajectory-B2	trajectory-A1	trajectory-A2	trajectory-A3
	trajectory-C1		6.57366228118406E-09	9.39120292642315E-09	0.041976035565175	0.00815491095821	0.015030842447119
	trajectory-B1	6.57366228118406E-09		0.793689131687568	7.1753212740111E-05	0.000724011746697	0.000523136902139
	trajectory-B2	9.39120292642315E-09	0.793689131687568		0.000111941306747	0.001154915701154	0.000885284265166
	trajectory-A1	0.041976035565175	7.1753212740111E-05	0.000111941306747		0.504157688012271	0.725345271990754
	trajectory-A2	0.00815491095821	0.000724011746697	0.001154915701154	0.504157688012271		0.530953268224535
	trajectory-A3	0.015030842447119	0.000523136902139	0.000885284265166	0.725345271990754	0.530953268224535	
	"""

	table = import_table_from_string(string, index = 'Trajectory')
	matrix = import_table_from_string(matrix_string, index = 'Index')
	print(matrix.to_string())
	distances = np.array([
		np.array([6.573662e-09, 9.391203e-09, 0.041976, 0.008155, 0.015031]),
		np.array([7.936891e-01, 0.000072, 0.000724, 0.000523]),
		np.array([0.000112, 0.001155, 0.000885]),
		np.array([0.504158, 0.725345]),
		np.array([0.530953])
	])
	#np.random.seed(4711)  # for repeatability of this tutorial
	#a = np.random.multivariate_normal([10, 0], [[3, 1], [1, 4]], size = [100, ])
	#b = np.random.multivariate_normal([0, 20], [[3, 1], [1, 4]], size = [50, ])
	#X = np.concatenate((a, b), )
	print(table.values)
	cpv = lambda a, b: (calculate_p_value(a, b, .03, .97).X)
	#cpv = partial(calculate_p_value, detected_cutoff = .03, fixed_cutoff = .97)
	Z = linkage(table.values, metric = cpv)
	print(Z)

	# calculate full dendrogram
	plt.figure(figsize = (10, 10))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')
	dendrogram(
		Z,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = table.index
	)
	#plt.show()
	plt.savefig('test_image.png')
