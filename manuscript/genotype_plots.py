import functools
import math
import random
from typing import Callable, Iterable, List, Optional, Union

import matplotlib.pyplot as plt
import numpy
import pandas

Number = Union[int, float]


def cdf(x, mean: float = 0.5, variance: float = 0.1, offset = 0):
	x += offset
	A = (x - mean) / (variance * math.sqrt(2))
	result = 0.5 * (1 + math.erf(A))

	return result


def degenerate_distribution(x: Number, point: Number, magnitude: Number) -> float:
	""" Will only return two values (zero and `magnitude` depending on whether the `x` values is less than `point`"""

	if x < point:
		return 0
	else:
		return magnitude


class Series:
	"""
	Parameters
	----------
	midpoint:float
		Determines where the midpoint of the series-specific feature (ex maximum delta, maximum value, etc.) is located on the x-axis.
	width: float
		Describes the general shape of the series (similar to the varaince).
	point_y: The y-location of the series-specific feature (inflection point, maximum, etc.).
	"""

	def __init__(self, midpoint: float = 0.5, width: float = 0.1, point_y = None, duration = None):
		"""
		Parameters
		----------
		midpoint: Number
			The timepoint where the series reaches its midpoint (usually the mean) value.
		width: Number
			Generally describes the shape of the series. Usually correspondins to the variance.
		point_y:Optional[Number]
			Resizes the distribution to fit between [0,`maximum`]
		duration:Optional[int]
			The number pf timepoints to generate
		"""
		self.midpoint = midpoint
		self.duration = duration
		self.point_y: Optional[Number] = point_y
		self.length = 100  # Total number of x-values to generate
		self.include_error = False

		self.x = numpy.linspace(0, 1, self.length)
		self.function = self.generate_series(self.midpoint, width)

	def plot(self, axes: Optional[plt.Axes] = None) -> plt.Axes:
		if axes is None:
			fig, axes = plt.subplots(figsize = (12, 8))

		y = self.y(self.x)
		if self.include_error:
			y = [(random.uniform(-0.01, 0.01)+i) for i in y]
		axes.plot(self.x, y)

		axes.set_ylim(0,1)
		axes.set_xlim(0,1)

		axes.set_xlabel("Duration of experiment", fontsize = 12)
		axes.set_ylabel("Genotype Frequency", fontsize = 12)
		axes.set_title("Genotypes observed during a population evolution experiment", fontsize = 18)

		axes.spines['right'].set_visible(False)
		axes.spines['top'].set_visible(False)
		return axes

	@staticmethod
	def generate_series(midpoint: float, width: float) -> Callable:
		raise NotImplementedError

	def y(self, xs: Iterable[Number]) -> List[float]:
		return [self.function(i) for i in xs]


class FixedSeries(Series):

	def generate_series(self, mean: float, width: float) -> Callable:
		cdffunction = functools.partial(cdf, mean = mean, variance = width)

		if self.point_y:
			# scale the series so that `self.maximum` corresponds to the fixed value
			return lambda s: cdffunction(s) * self.point_y
		else:
			return cdffunction


class FixedImmediateSeries(Series):
	def generate_series(self, midpoint: float, width: float) -> Callable:
		f = functools.partial(degenerate_distribution, point = midpoint, magnitude = 1)

		if self.point_y:
			return lambda s: f(s) * self.point_y
		else:
			return f


class ConstantSeries(Series):
	def __init__(self, midpoint: float = 0.5, width: float = 0.1, point_y = None, duration = None):
		super().__init__(midpoint, width, point_y, duration)
		# Make sure the maximum value (the mean value of the series) is set.
		assert self.point_y

	def generate_series(self, midpoint: float, width: float) -> Callable:
		"""
			The midpoint value does not matter for this function, since it is relatively constant.
			`width` describes the desired variance, `self.maximum` is the mean value of the series.
		"""

		def function(x: float):
			delta = random.gauss(width / 2, width / 2)
			return self.point_y + delta

		return function

	def y(self, xs: Iterable[Number]) -> List[float]:
		windowsize = 10
		ys = super().y(xs)
		# smooth the series
		s = pandas.Series(ys).rolling(windowsize).mean()
		# The first 4 timepoints will be NaN, so add some filler values
		values = s.tolist()
		values = values[windowsize-1:2*windowsize-2] + values[windowsize-1:]
		return values


class LostSeries(Series):
	""" Rises in frequency early in the experiment, but ultimately looses out when another genotype sweeps"""
	def __init__(self, midpoint: float = 0.5, width: float = 0.1, point_y = None, duration = None):
		super().__init__(midpoint, width, point_y, duration)
		# Make sure the maximum value (the mean value of the series) is set.
		assert self.point_y
	def generate_series(self, midpoint: float, width: float) -> Callable:
		"""
		Parameters
		----------
		mean:Number
			The x value where the series reaxes its maximum value.
		"""
		from scipy.stats import poisson
		scale = self.point_y / 0.2

		def function(x: Number) -> float:
			value = poisson.pmf(4, (x-(midpoint/2)) * 40) * scale
			if math.isnan(value):
				value = 0
			return value

		return function

class BeforeAfterFixed(Series):
	pass

class DetectedOnce(Series):
	# Remember to mention both fixed and unfixed variants
	def generate_series(self, midpoint, float, width:float = None):
		return lambda s: (1 if s == midpoint else 0)


class BeganFixed(Series):
	pass

if __name__ == "__main__":

	allseries = [
		FixedSeries(midpoint = 0.3, width = 0.08),
		FixedSeries(midpoint = 0.35, width = 0.1),
		FixedSeries(midpoint = 0.5, width = 0.3)
	]

	ax = None
	for series in allseries:
		ax = series.plot(ax)

	plt.show()
