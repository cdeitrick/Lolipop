
import matplotlib.pyplot as plt
import pandas
from loguru import logger

from muller import widgets


class DerivativePlot:
	def __init__(self):
		pass

	# noinspection PyMethodMayBeStatic
	def run(self, left: pandas.Series, right: pandas.Series):
		detected_left, detected_right = widgets.get_valid_points(left, right, 0.03, 0.97, inner = True)

		diff_left = detected_left.diff()
		diff_right = detected_right.diff()
		diff_series = list(zip(diff_left, diff_right))
		fig, ax = plt.subplots(figsize = (15, 15))
		ax.scatter(diff_left.values, diff_right.values)

		ax.set_xlabel("left derivative")
		ax.set_ylabel("right derivative")

		ax.set_xlim(-1, 1)
		ax.set_ylim(-1, 1)

		ax.axhline(0)
		ax.axvline(0)
		ax.plot([-1, 1], [-1, 1])
		ax.plot([-1, 1], [1, -1])
		correlated = sum([distance_correlated(i) for i in diff_series[1:]])
		anticorrelated = sum([distance_anticorrelated(i) for i in diff_series[1:]])

		logger.info(f"Correlation score: {correlated}")
		logger.info(f"Anticorrelation score: {anticorrelated}")

		plt.show()


import math


def distance_correlated(point) -> float:
	x, y = point
	d = abs(x - y) / math.sqrt(2)
	return d


def distance_anticorrelated(point):
	x, y = point
	d = abs(x + y) / math.sqrt(2)
	return d
