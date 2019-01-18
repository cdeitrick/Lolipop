import pandas
from dtw import dtw
import math
def normalize(series: pandas.Series)->pandas.Series:
	mean = series.mean()
	sigma = series.var()
	if math.isnan(sigma) or len(series) <2:
		normalized_series = None
	else:
		normalized_series = (series - mean)/sigma
	return normalized_series

def calculate_dtw(left:pandas.Series, right:pandas.Series, normal = True)->float:

	full = pandas.DataFrame([left, right]).T
	if left.name == right.name:
		reduced_df = full
	else:
		to_remove = list()
		for index, row in full.iterrows():
			valid_a = 0.03 < row[0] < 0.97
			valid_b = 0.03 < row[1] < 0.97
			if not valid_a and not valid_b:
				to_remove.append(index)
		reduced_df = full.drop(to_remove)

	reduced_df = full
	left_series = reduced_df.iloc[:,0]
	right_series= reduced_df.iloc[:,1]
	if normal:
		left_series = normalize(left_series)
		right_series = normalize(right_series)

	if left_series is None or right_series is None:
		distance,cost_matrix, acc_cost_matrix, path = 300, None, None, None
	else:
		distance, cost_matrix, acc_cost_matrix, path= dtw(left_series.values, right_series.values, dist = lambda x, y: abs(x - y))
	if math.isnan(distance):
		distance = 100
	if left.name == right.name:
		distance = 0


	return distance, cost_matrix, acc_cost_matrix, path

if __name__ == "__main__":
	from import_data import import_table_from_string
	trajectory_table = """
		Trajectory	0	1	3	4	6	7	9	10	12
		1	0	0	0	0	0	0	0	1	0
		2	0	0	0	0	0.052	0	0	0	0
		6	0	0	0	0	0.175	0	0	0	0
		7	0	0	0	0	0	0	0	0	0.054
		8	0	0	0	0	0	0	0	0	0.062
		16	0	0.058	0	0	0	0	0	0	0
		21	0	0	0	0	0	0	0	0.06	0
		23	0	0	0	0	0	0	0	0	0.162
		39	0	0	0	0	0	0	0	0.055	0
		43	0	0	0	0	0	0	0	0.349	0
		65	0	0	0	0	0.251	1	1	1	1
		71	0	1	1	1	1	1	1	1	1
	"""
	table = import_table_from_string(trajectory_table, index = 'Trajectory')
	from dtw import dtw
	l2_norm = lambda x, y: (x - y)**2
	l3_norm = lambda x,y: abs(x-y)

	for index, row in table.iterrows():
		distance, cost_matrix, acc_cost_matrix, path = calculate_dtw(table.loc['71'], row)
		print(index,"\t", distance)

	import matplotlib.pyplot as plt

	distance, cost_matrix, acc_cost_matrix, path = calculate_dtw(table.loc['71'], table.loc['65'])
	if acc_cost_matrix is not None:
		from pprint import pprint
		pprint(list(list(float(f"{j:.2f}") for j in i) for i in cost_matrix))
		plt.imshow(acc_cost_matrix.T, origin = 'lower', cmap = 'gray', interpolation = 'nearest')
		plt.plot(path[0], path[1], 'w')
		plt.show()