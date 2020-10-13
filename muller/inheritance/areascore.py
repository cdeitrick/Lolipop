import math
from typing import *

import pandas
from loguru import logger
from shapely import geometry
from shapely.errors import TopologicalError

Number = Union[float, int]

from . import polygon


DEBUG = False

def calculate_area(series:pandas.Series)->float:
	values = list()

	for previous, current in zip(series.values[:-1], series.values[1:]):
		if previous < current:
			lower = previous
			higher = current
		else:
			lower = current
			higher = previous

		area = lower + ((higher - lower)/ 2)
		values.append(area)
	return sum(values)



def area_of_series(series: pandas.Series) -> float:
	"""
		Calculates the area under the curve of a given series. Assumes that the series
		can be considered a sequence of polygons where vertices correspond to the
		x-values and y-values in the series.
	"""
	points = polygon.get_vertices(series)
	collection_of_polygons:List[List[Tuple[float,float]]] = polygon.isolate_polygons(points)
	area = 0
	for poly in collection_of_polygons:
		element_area = geometry.Polygon(poly).area
		area += element_area
	return area


def calculate_common_area(left: pandas.Series, right: pandas.Series) -> float:
	""" Calculates |X \cap Y|"""
	# Use the pandas method for now. It is probably slower, but shouldn't add too much time, and should be more readable.

	left_points = polygon.get_vertices(left) if isinstance(left, pandas.Series) else left  # Only decompose if the input was not a list of points.
	right_points = polygon.get_vertices(right) if isinstance(right, pandas.Series) else right

	# Test if one of the series has regions where it is undetected between other areas where it is detected.
	# Basically whether the series needs to be represented by two separate points
	# Since polygon.separate() can handle the case where the series does not need to be modified, may as well run it on all
	# series instead of testing for disjoint polygons first.
	left_collection = polygon.isolate_polygons(left_points)
	right_collection = polygon.isolate_polygons(right_points)

	left_collection_polygon = [geometry.Polygon(i) for i in left_collection if len(i) > 2]
	right_collection_polygon = [geometry.Polygon(i) for i in right_collection if len(i) > 2]

	left_collection_polygon = geometry.MultiPolygon(left_collection_polygon)
	right_collection_polygon = geometry.MultiPolygon(right_collection_polygon)
	try:
		intersection = left_collection_polygon.intersection(right_collection_polygon)
	except TopologicalError as exception:
		logger.error(f"Cannot calculate the intersection:")
		logger.error(f"\t{left_collection_polygon}")
		logger.error(f"\t{right_collection_polygon}")
		raise exception

	total_intersection_area = intersection.area

	return total_intersection_area


def X_and_Y_polygon(left: geometry.MultiPolygon, right: geometry.MultiPolygon) -> float:
	try:
		return left.intersection(right).area
	except TopologicalError:
		logger.error(f"Polygonal error")
		return 0


def X_or_Y_numeric(left: float, right: float, x_and_y: float) -> float:
	""" Calculates the OR area between `left` and `right`
		Parameters
		----------
		left, right: pandas.Series
			The two series to calculate the exclusive area on.
		x_and_y: float
			Precomputed values for the area of `left`, `right`, and x_and_y.
	"""
	# Need to subtract `x_and_y` since it is counted twice in `left` and `right`
	return left + right - x_and_y


def X_or_Y(left: pandas.Series, right: pandas.Series) -> float:
	""" Calculates the OR area between `left` and `right`
		Parameters
		----------
		left, right: pandas.Series
			The two series to calculate the exclusive area on.
	"""

	left_polygons = polygon.as_polygon(left)
	right_polygons = polygon.as_polygon(right)

	multipolygon = geometry.MultiPolygon([left_polygons, right_polygons])
	area = multipolygon.area - left_polygons.intersection(right_polygons).area

	return area


def X_xor_Y(left: pandas.Series, right: pandas.Series) -> float:
	left_polygons = polygon.as_polygon(left)
	right_polygons = polygon.as_polygon(right)

	xor_polygons = left_polygons.symmetric_difference(right_polygons)
	return xor_polygons.area


def difference_polygon(left: geometry.MultiPolygon, right: geometry.MultiPolygon) -> float:
	""" Returns the area of `left` not in `right`"""
	try:
		area = left.difference(right).area
	except TopologicalError:
		logger.error(f"Polygonal error")
		area = 0
	return area


def jaccard_distance_numeric(x_or_y: float, x_and_y: float) -> float:
	""" Calculated the jaccard distance between two series using pre-computed area values.
		Parameters
		----------
		x_or_y, x_and_y: float
			precomputed x_or_y and x_and_y, respectively.
	"""
	return (x_or_y - x_and_y) / x_or_y


def jaccard_distance(left: pandas.Series, right: pandas.Series) -> float:
	""" Calulated the jaccard distance between two series.
		Parameters
		----------
		left, right: pandas.Series
			The two genotypes to calculate the jaccard distance on.
	"""

	x_or_y = X_or_Y(left, right)
	x_and_y = X_and_Y(left, right)
	return jaccard_distance_numeric(x_or_y, x_and_y)


def jaccard_subset_numeric(left: float, right: float) -> float:
	""" Operates on pre-computed areas to help reduce redundant calculations.
		Calculates the jaccard distance assuming `right` is a subset of `left`.
	Parameters
	----------
	left, right: pandas.Series
	left, right: float
	"""
	return (left - right) / left


def jaccard_subset(left: pandas.Series, right: pandas.Series) -> float:
	"""
		Calculates the jaccard distance assuming `right` is a subset of `left`.
	Parameters
	----------
	left, right: pandas.Series
	left, right: float
	"""
	left_area = area_of_series(left)
	right_area = area_of_series(right)
	return jaccard_subset_numeric(left_area, right_area)


def is_subset(left: pandas.Series, right: pandas.Series) -> bool:
	""" A rewrite that uses the shapely library instead."""
	left_poly = polygon.as_polygon(left)
	right_poly = polygon.as_polygon(right)

	try:
		area_intersection = left_poly.intersection(right_poly).area
	except TopologicalError:
		logger.error(f"{left.values}, {left.index}")
		logger.error(f"{right.values}, {right.index}")
		logger.error(f"Left polygon: {left_poly}")
		logger.error(f"Right polygon: {right_poly}")
		logger.error(f"TopologyError")
		return False

	area_right = right_poly.area

	result = math.isclose(area_intersection, area_right, abs_tol = 0.03 ** 2)  # Tolerance is square of the detection limit.
	return result


def is_subset_polygon(left: geometry.MultiPolygon, right: geometry.MultiPolygon) -> bool:
	area_union = left.union(right).area
	area_intersection = left.intersection(right).area
	area_left = left.area
	area_right = right.area

	result = math.isclose(area_intersection, area_right, abs_tol = 0.03 ** 2)

	jaccard_expected = (area_left - area_right) / area_left
	jaccard_actual = (area_union - area_intersection) / area_union

	if DEBUG:
		logger.debug(f"left, right -> {area_left}, {area_right}")
		logger.debug(f"Je = ({area_left:.2f} - {area_right:.2f}) / {area_left:.2f} = {jaccard_expected:.2f}")
		logger.debug(f"Ja = ({area_union:.2f} - {area_intersection:.2f}) / {area_union:.2f} = {jaccard_actual:.2f}")

		logger.debug(f"{jaccard_expected} == {jaccard_actual} -> {result}")
	result = math.isclose(jaccard_expected, jaccard_actual, abs_tol = 0.1)
	return result


X_and_Y = calculate_common_area



