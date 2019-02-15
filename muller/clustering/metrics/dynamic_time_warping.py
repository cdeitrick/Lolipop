import math

import pandas
from dtw import dtw


def normalize(series: pandas.Series) -> pandas.Series:
	mean = series.mean()
	sigma = series.var()
	if math.isnan(sigma) or len(series) < 2:
		normalized_series = None
	else:
		normalized_series = (series - mean) / sigma
	return normalized_series


def filter_invalid_timepoints(df: pandas.DataFrame) -> pandas.DataFrame:
	to_remove = list()
	for index, row in df.iterrows():
		valid_a = 0.03 < row[0] < 0.97
		valid_b = 0.03 < row[1] < 0.97
		if not valid_a and not valid_b:
			to_remove.append(index)
	return df.drop(to_remove)


def calculate_dtw(left: pandas.Series, right: pandas.Series, normal = False) -> float:
	full = pandas.DataFrame([left, right]).T
	reduced_df = full
	left_series = reduced_df.iloc[:, 0]
	right_series = reduced_df.iloc[:, 1]
	if normal:
		left_series = normalize(left_series)
		right_series = normalize(right_series)
	if left_series is None or right_series is None:
		distance, cost_matrix, acc_cost_matrix, path = 100, None, None, None
	else:
		mean_series = reduced_df.mean(axis = 1)
		sigma = mean_series.var()
		distance_metric = lambda x, y: abs(x - y) / sigma

		distance, cost_matrix, acc_cost_matrix, path = dtw(left_series.values, right_series.values, dist = distance_metric, w = 3, s = 2)

	if left.name == right.name:
		distance = 0

	# return distance, cost_matrix, acc_cost_matrix, path
	return distance


if __name__ == "__main__":
	pass
