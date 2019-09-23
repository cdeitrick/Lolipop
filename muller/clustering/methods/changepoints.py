import math
from typing import List

import matplotlib.pyplot as plt
import pandas
from scipy import stats


def detect_changepoint_cbs(values: List[int]):
	alpha = 0.05
	statistics = list()
	for index in range(0, len(values)):
		segment = values[:index]
		Zi = calculate_Zi(len(segment), sum(segment), len(values), sum(values))
		statistics.append((index - 1, abs(Zi)))
	Zb = max(statistics, key = lambda s: s[1])
	return Zb, statistics


def calculate_Zi(i, Si, n, Sn):
	A = (Si / i) - ((Sn - Si) / (n - i))
	B = (1 / i) + (1 / (n - i))

	return A / math.sqrt(B)


def get_maximums(values):
	series = pandas.Series(values)
	diff_series = series.diff()
	rank = diff_series.rank(method = 'first')
	maximum_difference = diff_series.max()
	diff_series = diff_series[diff_series > (maximum_difference)]


def detect_changepoint(values: List[float]):
	# remove 0's
	values = [i for i in values if i]
	# Remove the maximum values, which were used for overlapping series
	values = [i for i in values if not math.isclose(i, max(values))]

	statistics = list()
	for index, value in enumerate(values):
		if index < 3: statistics.append((index, 0))
		else:
			segment = values[:index]
			point = values[index]
			statistic, pvalue = stats.ttest_1samp(segment, point)
			statistics.append((index, abs(statistic)))

	index, statistic = max(statistics, key = lambda s: s[1])
	return (index, values[index]), statistics


def table_to_list(table: pandas.DataFrame):
	values = list()
	for index, row in table.iterrows():
		values += list(row.values)
	values = sorted(values)
	# Remove 0 distance values
	values = [i for i in values if i]
	# Remove the maximum values, which were used for overlapping series
	values = [i for i in values if not math.isclose(i, max(values))]
	return values


def plot_normality(values: List[float]):
	import seaborn
	seaborn.distplot(values, kde = False, color = "b")
	plt.show()


def plot_changepoints(values: List[float], changepoint: int, statistics = None):
	fig, ax = plt.subplots(figsize = (20, 10))
	x = list(range(len(values)))
	ax.scatter(x, values)
	ax.axvline(changepoint)
	if statistics is not None:
		x, y = zip(*statistics)
		ax.scatter(x, y)
	series = pandas.Series(values, index = x)
	series = series.diff()
	ax.scatter(series.index, series.values, c = 'red')
	plt.show()
