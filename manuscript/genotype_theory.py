import math
from typing import Iterable, List, Union, Any

import numpy
import pandas
from loguru import logger

Number = Union[int, float]
from pathlib import Path

def sumlist(items):
	r = list()
	for i in items:
		r += i
	return r
# Try to design a system which adds `features` of series together rather than trying to generate series as a whole.
def cdf(x, mean: float = 0.5, variance: float = 0.1, offset = 0):
	x += offset
	A = (x - mean) / (variance * math.sqrt(2))
	result = 0.5 * (1 + math.erf(A))

	return result


def newrange(values, minimum: Number, maximum: Number):
	obs_min = min(values)
	obs_max = max(values)
	obs_delta = obs_max - obs_min
	delta = maximum - minimum

	convert = lambda s: ((s - obs_min) * delta / obs_delta) + minimum

	return [convert(i) for i in values]


class Base:
	def __init__(self, *args, **kwargs):
		self.value = []
		self.length = 10

	def __str__(self) -> str:
		string = f"{self.__name__}({self.value})"
		return string

	def __add__(self, other: 'Base') -> 'Values':
		return Values(other.get() + self.get())

	def get(self, length: int) -> List[Number]:
		pass


class Values(Base):
	def __init__(self, value: Iterable[Number]):
		super().__init__()
		self.value = value

	def get(self, length: int = None):
		return self.value


class Constant(Base):
	def __init__(self, value: float):
		super().__init__()
		self.value = value

	def get(self, length: int) -> List[Number]:
		return [self.value for _ in range(length)]


class Inflection(Base):
	""" Generates an inflection series from `min` to `max`"""

	def __init__(self, minimum, maximum):
		super().__init__()
		self.minimum = minimum
		self.maximum = maximum

	def get(self, length: int):
		domain = numpy.linspace(-3, 3, length)
		values = [0.5 * (1 + math.erf(i)) for i in domain]

		coerced_values = newrange(values, self.minimum, self.maximum)

		return coerced_values

	def get_series(self, length: int, position: float):
		inflection_length = int(length * 0.2)
		remaining_length = length - inflection_length
		start_length = int(remaining_length * position)
		end_length = remaining_length - start_length

		return Constant(self.minimum).get(start_length) + self.get(inflection_length) + Constant(self.maximum).get(end_length)


class Hill(Base):
	def __init__(self, height: float):
		super().__init__()
		self.height = height
		self.inflection_width = 0.1

	def get(self, length: int) -> List[Number]:
		start = int(length / 2)
		end = length - start
		left = Inflection(0, self.height).get(start)
		right = Inflection(0, self.height).get(end)[::-1]
		hill_segment =  left + right

		return hill_segment

	def get_series(self, length: int, position: float, width: float):
		hill_length = int(length * width)
		hill = self.get(hill_length)

		zero = Constant(0)
		remaining_length = length - hill_length
		start = int(remaining_length * position)
		end = remaining_length - start

		return zero.get(start) + hill + zero.get(end)


def totable(series: List[List[Number]], labels: List[str] = None):
	if not labels:
		labels = range(1, len(series)+1)
	labels = [f'genotype-{i}' for i in labels]
	columns = numpy.linspace(1, 100, len(series[0]))
	table = pandas.DataFrame(series, index = labels, columns = [str(int(i)) for i in columns])
	table.index.name = "Genotype"
	return table


def modelplot(table: pandas.DataFrame, filename: Union[str, Path] = None, colors:List[str] = None):
	import matplotlib.pyplot as plt
	fig, ax = plt.subplots(figsize = (12, 8))

	for index, (label, row) in enumerate(table.iterrows()):
		color = colors[index%len(colors)] if colors else None

		ax.plot(row.index, row.values, c = color, label = label, linewidth = 12)
	if filename:
		plt.savefig(filename)
	else:
		plt.show()


def model_periodic_selection():
	inflection = Inflection(0, 1)
	zero = Constant(0)
	one = Constant(1)
	results = [
		zero.get(3) + inflection.get(22) + one.get(75),
		zero.get(25) + inflection.get(20) + one.get(55),
		zero.get(30) + inflection.get(20) + one.get(50),
		zero.get(50) + inflection.get(50)
	]
	table = totable(results)
	table.to_csv("periodic_selection.tsv")
	modelplot(table, "periodic_selection.png")


def model_clonal_interferance():
	total_length = 100
	inflection = Inflection(0, 1)
	hill = Hill(0.3)
	orange_segment = [
		Constant(0).get(20),
		Inflection(0, 0.6).get(20),
		Constant(0.6).get(10),
		Inflection(0.6, 0.7).get(10),
		Constant(0.7).get(15),
		Inflection(0.7, 1).get(10),
		Constant(1).get(5)
	]
	blue_segment = [
		Constant(0).get(21),

		Inflection(0, 0.4).get(18),
		Constant(0.4).get(11),
		Inflection(0.3, 0.4).get(10)[::-1],

		Constant(0.3).get(10),
		Inflection(0, 0.3).get(10)[::-1],
		Constant(0).get(20)
	]

	purple_segment = [
		Constant(0).get(50),
		Inflection(0, 0.7).get(15),
		Constant(0.7).get(10),
		Inflection(0.7, 1).get(10),
		Constant(1).get(15)
	]
	green_segment = [
		Constant(0).get(50),
		Inflection(0,0.3).get(10),
		Constant(0.3).get(10),
		Inflection(0,0.3).get(10)[::-1],
		Constant(0).get(20)
	]

	results = [
		inflection.get_series(total_length, 0.1),
		#hill.get_series(total_length, 0.7, 0.15),
		sumlist(orange_segment),
		sumlist(blue_segment),
		sumlist(purple_segment),
		sumlist(green_segment),
		Inflection(0,1).get_series(100, 0.85)
	]
	colors = ['mauve', 'orange', 'blue', 'purple', 'green', 'red']
	table = totable(results, labels = colors)
	table.to_csv("clonal_interferance.tsv", sep = "\t")
	colors[0] = 'purple'
	modelplot(table, "clonal_interferance.png", colors = colors)

def model_strong_selection():
	total_length = 100

	orange_segment = [
		Constant(0).get(10),
		Inflection(0,.2).get(20),
		Constant(0).get(70)
	]

	blue_segment = [
		Constant(0).get(27),
		Inflection(0,1).get(10),
		Constant(1).get(63)
	]

	green_segment = [
		Constant(0).get(50),
		Inflection(0,.5).get(25),
		Constant(.5).get(5),
		Constant(1).get(20)
	]

	purple_segment = [
		Constant(0).get(50),
		Inflection(0,.4).get(30),
		Constant(0).get(20)
	]

	red_segment = [
		Constant(0).get(80),
		Inflection(0,1).get(5),
		Constant(1).get(15)
	]

	results = [
		Inflection(0,1).get_series(total_length, 0.1),
		sumlist(orange_segment),
		sumlist(blue_segment),
		sumlist(green_segment),
		sumlist(purple_segment),
		sumlist(red_segment)
	]
	colors = ['mauve', 'orange', 'blue', 'green', 'purple', 'red']
	table = totable(results, labels = colors)
	table.to_csv('model_strong_selection.tsv', sep = "\t")
	colors[0] = 'purple'
	modelplot(table, "model_strong_selection.png",colors = colors)
if __name__ == "__main__":
	model_clonal_interferance()

