import math

import pandas
from loguru import logger

try:
	from muller.inheritance.areascore import area_of_series, calculate_common_area
except ModuleNotFoundError:
	try:
		from ...inheritance.areascore import area_of_series, calculate_common_area
	except ValueError:
		from muller.inheritance.areascore import area_of_series, calculate_common_area


def minkowski_distance(left: pandas.Series, right: pandas.Series, p: int = 2) -> float:
	""" Calculates the minkowski distance between two series. Essentially just a generic lp-norm.
		Parameters
		----------
		left:pandas.Series
		right:pandas.Series
		p:int
	"""
	total = sum([math.pow(abs(i - j), p) for i, j in zip(left.tolist(), right.tolist())])
	return math.pow(total, 1 / p)


def pearson_correlation_distance(left: pandas.Series, right: pandas.Series, adjusted = True) -> float:
	"""
		Calculates the pearson correlation between two series. The resulting value lies in the range [0,2].
	Parameters
	----------
	left:pandas.Series
	right:pandas.Series
	adjusted: bool; default True
		Whther to adjust the coefficient for small sample sizes.

	Returns
	-------
	float
	"""

	pcc = left.corr(right, method = 'pearson')
	# Adjust due to sample size
	if adjusted:
		adjusted_pcc = adjust_correlation_coefficient(pcc, len(left))
	else:
		adjusted_pcc = pcc
	# convert to distance metric.
	return 1 - adjusted_pcc


def adjust_correlation_coefficient(r: float, n: int) -> float:
	value = (1 - r ** 2) / (2 * n)
	ra = r * (1 + value)
	return ra


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
	n = len(mean)
	sigma_freq: pandas.Series = mean.mul(1 - mean)
	# Difference of frequencies at each timepoint
	# difference: pandas.Series = not_detected_fixed_df.iloc[:, 0] - not_detected_fixed_df.iloc[:, 1]
	difference = not_detected_fixed_df.diff(axis = 1).iloc[:, 1]

	sigma_pair: float = sigma_freq.sum() / n ** 2
	# Sum of differences
	difference_mean: float = abs(difference).sum() / n

	X = difference_mean / (math.sqrt(2 * sigma_pair))

	return X


def binomial_probability(left: pandas.Series, right: pandas.Series) -> float:
	""" Calculates the probability that the distance between `left` and `right` can be explained by experimental error.
		Note that higher probabilities indicate a closer relationship between two series in contrast to distance measurements where
		lower distances correspond to better matches. This can be accounted for by subtracting the probability value from 1.

	"""
	X = binomial_distance(left, right)
	value = 1 - math.erf(X)
	logger.error("Set a patch so that binomialp represents probability")
	value = 1- value
	return value


def jaccard_distance(left: pandas.Series, right: pandas.Series) -> float:
	area_left = area_of_series(left)
	area_right = area_of_series(right)
	area_shared = calculate_common_area(left, right)
	j = area_shared / (area_left + area_right - area_shared)
	return 1 - j


def calculate_distance(left: pandas.Series, right: pandas.Series, metric: str) -> float:
	if metric == 'binomial':
		# default option so placed first to minimize equality calculations.
		distance_between_series = binomial_distance(left, right)
	elif metric == 'pearson':
		distance_between_series = pearson_correlation_distance(left, right)
	elif metric == 'minkowski':
		distance_between_series = minkowski_distance(left, right)
	elif metric == 'similarity' or metric == 'binomialp':
		distance_between_series = binomial_probability(left, right)
	elif metric == 'jaccard':
		distance_between_series = jaccard_distance(left, right)
	elif metric == "combined":
		distance_between_series_pearson = pearson_correlation_distance(left, right)
		distance_between_series_minkowski = minkowski_distance(left, right, 2)
		distance_between_series = (2 * distance_between_series_pearson) + distance_between_series_minkowski
	else:
		message = f"'{metric}' is not an available metric."
		raise ValueError(message)
	return distance_between_series
