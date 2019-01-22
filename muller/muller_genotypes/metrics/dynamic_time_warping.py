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

def filter_invalid_timepoints(df:pandas.DataFrame)->pandas.DataFrame:
	to_remove = list()
	for index, row in df.iterrows():
		valid_a = 0.03 < row[0] < 0.97
		valid_b = 0.03 < row[1] < 0.97
		if not valid_a and not valid_b:
			to_remove.append(index)
	return df.drop(to_remove)

def calculate_dtw(left:pandas.Series, right:pandas.Series, normal = False)->float:

	full = pandas.DataFrame([left, right]).T
	if left.name == right.name:
		reduced_df = full
	else:
		reduced_df = filter_invalid_timepoints(full)
	reduced_df = full
	left_series = reduced_df.iloc[:,0]
	right_series= reduced_df.iloc[:,1]
	if normal:
		left_series = normalize(left_series)
		right_series = normalize(right_series)
	if left_series is None or right_series is None:
		distance, cost_matrix, acc_cost_matrix, path = 100, None, None, None
	else:
		mean_series = reduced_df.mean(axis = 1)
		sigma = mean_series.var()
		distance_metric = lambda x, y: abs(x - y)/sigma

		distance, cost_matrix, acc_cost_matrix, path= dtw(left_series.values, right_series.values, dist = distance_metric, w = 3, s = 2)

	if left.name == right.name:
		distance = 0


	return distance, cost_matrix, acc_cost_matrix, path

if __name__ == "__main__":
	from import_data import import_table_from_string
	from muller_genotypes.metrics.similarity import calculate_p_value
	trajectory_table = """
		Trajectory	0	17	25	44	66	75	90
		1	0	0	0.261	1	1	1	1
		20	0	0	0	0.138	0.295	0	0.081
		4	0	0	0	0	0.211	0.811	0.813
		8	0	0	0	0	0.345	0.833	0.793
		16	0	0	0	0	0.209	0.209	0
		13	0	0	0	0	0.258	0.057	0.075
		15	0	0	0.066	0.104	0.062	0	0
		11	0	0	0	0.108	0.151	0	0
		17	0	0	0	0	0	0.266	0.312
		9	0	0	0	0	0	0.269	0.34
		18	0	0	0	0.115	0	0.131	0
		21	0	0	0	0.114	0	0.11	0.123
		14	0	0.38	0.432	0	0	0	0
		6	0	0	0	0	0	1	1
		2	0	0	0	0.525	0.454	0.911	0.91
		3	0	0	0	0.147	0.45	0.924	0.887
		7	0	0	0	0.273	0.781	1	1
		19	0	0	0	0.188	0.171	0.232	0.244
		5	0	0	0	0.403	0.489	0.057	0.08
		10	0	0	0.117	0	0	0	0.103
	"""
	table = import_table_from_string(trajectory_table, index = 'Trajectory')
	from dtw import dtw
	l2_norm = lambda x, y: (x - y)**2
	l3_norm = lambda x,y: abs(x-y)

	for index, row in table.iterrows():
		distance, cost_matrix, acc_cost_matrix, path = calculate_dtw(table.loc['1'], row, False)
		p_value = calculate_p_value(table.loc['1'], row, .03, .97).X
		print(f"{index}\t{distance:.2f}\t{p_value:.2f}")

	import matplotlib.pyplot as plt
	red = table.loc['1']
	blue = table.loc['6']

	if acc_cost_matrix is not None:
		from pprint import pprint
		#pprint(list(list(float(f"{j:.2f}") for j in i) for i in cost_matrix))
		plt.imshow(acc_cost_matrix.T, origin = 'lower', cmap = 'gray', interpolation = 'nearest')
		plt.plot(path[0], path[1], 'w')
		plt.show()

	redn = normalize(red)
	bluen = normalize(blue)
	distance, cost_matrix, acc_cost_matrix, path = calculate_dtw(red,blue)
	print(distance)
	plt.plot(red.index, red.values, color = 'r')
	plt.plot(redn.index, redn.values, color = 'r')
	plt.plot(blue.index, blue.values, color = 'b')
	plt.plot(bluen.index, bluen.values, color = 'b')
	plt.show()


