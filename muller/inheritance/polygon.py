from typing import Any, Iterable, List, Tuple, Union

import pandas
from shapely import geometry

Number = Union[float, int]
PointType = Tuple[Number, Number]
import collections
import itertools


def sliding_window_legacy(iterable: Iterable, width: int):
	for index in range(len(iterable)):
		try:
			if index + width > len(iterable):
				raise StopIteration  # Because otherwise will still return last last n elements.
			yield tuple(iterable[index:index + width])
		except (IndexError, StopIteration):
			break


def sliding_window(seq, n):
	""" A sequence of overlapping subsequences
	>>> list(sliding_window(2, [1, 2, 3, 4]))
	[(1, 2), (2, 3), (3, 4)]
	This function creates a sliding window suitable for transformations like
	sliding means / smoothing
	>>> mean = lambda seq: float(sum(seq)) / len(seq)
	>>> list(map(mean, sliding_window(2, [1, 2, 3, 4])))
	[1.5, 2.5, 3.5]
	"""
	return zip(*(collections.deque(itertools.islice(it, i), 0) or it
		for i, it in enumerate(itertools.tee(seq, n))))


def get_neighbor(index: pandas.Index, item: Any, previous: bool = True) -> Any:
	"""
		Return the neighboring value to `item`.
	Parameters
	----------
	index: The full index of the series. NOT just the index of detected values
	item
	previous

	Returns
	-------

	"""
	# Get the integer position of the item
	index = list(index)  # pandas.Index uses circular indexing, so the last key is treated as being prior to the initial key.
	offset = -1 if previous else 1

	location = index.index(item)
	if location != 0:
		location += offset
	try:
		return index[location]
	except IndexError:
		return item


def _decompose_correct_split_series(series: pandas.Series):
	""" Corrects the case where a series is undetected at a single point, but otherwise is detected.
		This will cause a topology error in shapely.
	exmple:
		[0.97, 0.0, 0.97, 0.97, 0.97, 0.87, 0.97]
	"""
	default_value = 0.001  # Should be small enough to not affect much.
	# Figure out which points need to be modified.
	indicies = list()

	for index, window in enumerate(sliding_window(series, 3)):
		if window[0] > 0 and window[2] > 0 and window[1] == 0:
			indicies.append(index)
	for i in indicies:
		# guaranteed to not be the last index due to how sliding_window() works.
		series[i + 1] = default_value
	return series


def _correct_one_dimensional_sections(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
	""" Removes adjacent 0's representing a one-dimensional offshoot form the polygon."""

	for index, (x, y) in enumerate(points):
		if y == 0:
			try:
				xp, yp = points[index - 1]
				if yp == 0:
					points[index - 1] = (xp, 0.001)
			except IndexError:
				pass
	return points



def decompose(series: pandas.Series, use_index:bool = False) -> List[PointType]:

	# Make sure one-dimensional points are omitted.
	minimum = 0.0001
	masked_series = series.mask(lambda s: s < minimum, minimum).tolist()
	masked_series[0] = series.values[0]
	masked_series[-1] = series.values[-1]

	if use_index:
		x_values = series.index
	else:
		x_values = list(range(len(series)))

	points = list(zip(x_values, masked_series))
	if points: # Check if the list is empty
		first_point = points[0]
		if first_point[1] != 0:
			points = [(first_point[0], 0)] + points

		last_point = points[-1]
		if last_point[1] != 0:
			points.append((last_point[0], 0))
	return points


def separate_main(series: List[PointType]) -> List[List[PointType]]:
	""" Separates a series into segregated polygons assuming there is a region of no detection between them."""
	if not is_multiple(series):
		return [series]
	for index, element in enumerate(sliding_window(series, 2), start = 1):  # Start indexing at 1 since starting with second element.
		if element[0][1] == 0 and element[1][1] == 0:
			# Separate the polygons to the left and right of this region.
			left = series[:index]
			right = series[index:]
			return separate_main(left) + separate_main(right)


def separate(series: List[PointType]) -> List[List[PointType]]:
	""" Implemented here until I refactor the above function to omit small sequences"""
	return [i for i in separate_main(series) if len(i) > 2]


def is_multiple(series: List[Tuple[Number, Number]]):
	# Test if one of the series has regions where it is undetected between other areas where it is detected.
	# Basically whether the series needs to be represented by two separate points
	for window in sliding_window(series, 2):
		if window[0][1] == 0 and window[1][1] == 0:
			return True
	return False


def as_polygon(series: pandas.Series) -> geometry.MultiPolygon:
	""" Converts a series object to a polygon object.
		[0.97 0.   0.97 0.97 0.97 0.87 0.97] Fails
	"""
	points = decompose(series)
	result = geometry.Polygon(points)
	return result
