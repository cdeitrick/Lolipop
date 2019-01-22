import math

import pandas
from dtw import dtw


def minkowski_distance(left: pandas.Series, right: pandas.Series, p: int) -> float:
	""" Calculates the minkowski distance between two series. Essentially just a generic lp-norm.
		Parameters
		----------
		left:pandas.Series
		right:pandas.Series
		p:int
	"""
	total = sum([math.pow(abs(i - j), p) for i, j in zip(left.tolist(), right.tolist())])
	return math.pow(total, 1 / p)


def pearson_correlation_distance(left: pandas.Series, right: pandas.Series) -> float:
	"""
		Calculates the pearson correlation between two series. The resulting value lies in the range [0,2].
	Parameters
	----------
	left:pandas.Series
	right:pandas.Series

	Returns
	-------
	float
	"""
	covariance = left.cov(right)
	den = left.std() * right.std()
	# The pearson correlation coefficient is in the range [-1, 1], where values closer to 1 indicate perfect similarity
	pcc = covariance / den
	# convert to distance metric.
	return 1 - pcc


def dynamic_time_warping(left: pandas.Series, right: pandas.Series) -> float:
	""" Calculates the DTW distance between two series.
	"""
	cost_metric = lambda x, y: (x - y) ** 2
	distance, cost_matrix, acc_cost_matrix, path = dtw(left.values, right.values, dist = cost_metric, w = 3, s = 2)
	return distance


# noinspection PyTypeChecker
def binomial_distance(left: pandas.Series, right: pandas.Series) -> float:
	""" Based on the binomial calculations present in the original matlab scripts."""
	# Find the mean frequency of each timepoint
	# index is timepoints,  values are frequencies
	not_detected_fixed_df = pandas.concat([left, right], axis = 1)
	mean: pandas.Series = not_detected_fixed_df.mean(axis = 1)

	# Calculate sigma_freq
	# E(sigma) = (1/n) sum(sigma) = (1/n) sum(np(1-p)) == sum(p(1-p)
	# E(sigma_p) = (1/n) E(sigma) == 1/n(sum(p(1-p))
	# E(d_bar) = 1/n(sum(di)) == 1/n (n*sum(di))
	# pandas.Series.radd is slow for some reason. Use '-' operator instead.
	sigma_freq: pandas.Series = mean.mul(1 - mean)
	# Difference of frequencies at each timepoint
	# difference: pandas.Series = not_detected_fixed_df.iloc[:, 0] - not_detected_fixed_df.iloc[:, 1]
	difference = not_detected_fixed_df.diff(axis = 1).iloc[:, 1]
	sigma_pair: float = sigma_freq.sum() / len(mean)
	# Sum of differences
	difference_mean: float = abs(difference).sum()

	X = difference_mean / (math.sqrt(2 * sigma_pair))

	return X


def binomial_probability(left: pandas.Series, right: pandas.Series) -> float:
	X = binomial_distance(left, right)
	value = 1 - math.erf(X)
	return value


if __name__ == "__main__":
	from import_data import import_table_from_string

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

	left = table.loc['7']
	right = table.loc['4']
	center = table.loc['3']
	print(minkowski_distance(left, right, 2))
	print(minkowski_distance(left, center, 2))

	print("Pearson")
	print(pearson_correlation_distance(left, right))
	print(pearson_correlation_distance(left, center))

	print("DTW")
	print(dynamic_time_warping(left, right))
	print(dynamic_time_warping(left, center))
