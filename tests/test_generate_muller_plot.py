import pytest
from muller.graphics.generate_muller_plot import *
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