import pytest

from muller import dataio
from muller.graphics.muller_plot import *


def test_relocate_point_simple():
	point = (.1, .613)
	locations = [(.1, .626), (5, .35)]

	expected_point = (.6, .513)
	result_point = relocate_point(point, locations)
	assert pytest.approx(expected_point, result_point)


def test_relocate_point_complex():
	point = (0.1, 0.426)
	locations = [
		(1.0, 0.154),
		(6.0, 0.0821442859047873),
		(0.1, 0.3475000000000004),
		(1.0, 0.26),
		(1.0, 0.4985),
		(4.5, 0.8981050818260126),
		(0.1, 0.6129285714285719)
	]
	expected_point = (.6, .726)
	result_point = relocate_point(point, locations)

	assert pytest.approx(expected_point, result_point)


def test_get_font_properties():
	expected_result = {
		'size':  9,
		'color': '#333333'
	}

	assert get_font_properties("#99d8c9") == expected_result

	assert get_font_properties("#00441b")['color'] == '#FFFFFF'


def test_find_closest_point():
	point = (10, 15)

	points = [(0, 1), (10, 10), (11, 11), (12, 12)]
	assert find_closest_point(point, points) == (12, 12)


def test_unique_everseen():
	assert list(unique_everseen("AABBCCCADAACCFD")) == "A B C D F".split()


def test_generate_muller_series():
	string = """
	Generation	Group_id	Frequency
	0	genotype-1	0
	0	genotype-12	0
	2	genotype-1	.1
	2	genotype-12	1
	4	genotype-1	.4
	4	genotype-12	.9
	6	genotype-1	1
	6	genotype-12	.5
	"""
	table = dataio.import_table(string)
	palette = {'genotype-1': '#AAAAAA', 'genotype-12': '#BBBBBB'}

	expected_x = [0, 2, 4, 6]
	expected_y = [
		[0, .1, .4, 1], [0, 1, .9, .5]
	]
	expected_colors = ['#AAAAAA', '#BBBBBB']
	expected_labels = ['genotype-1', 'genotype-12']

	result_x, result_y, result_colors, result_labels = generate_muller_series(table, palette)

	assert result_x == expected_x
	assert result_y == expected_y
	assert result_colors == expected_colors
	assert result_labels == expected_labels
