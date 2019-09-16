from pathlib import Path

import pandas
import pytest

from muller import dataio
from muller.inheritance import polygon, areascore
FOLDER_DATA = Path(__file__).parent.parent / "data"


@pytest.fixture
def trajectory_table() -> pandas.DataFrame:
	filename_table = FOLDER_DATA / "truthsets" / "truthset.model.area.xlsx"
	return dataio.import_table(filename_table, sheet_name = "data", index = 'Trajectory')


@pytest.mark.parametrize(
	"key, expected",
	[
		('A', [(0, 0), (1, .1), (2, .1), (3, .2), (4, .2), (5, .3), (6, .3), (6, 0)]),
		('B', [(0, 0)] + [(i, .1) for i in range(7)] + [(6, 0)]),
		('C', [(0, 0), (0, 1), (1, .9), (2, .8), (3, .7), (4, .6), (5, .5), (6, .4), (6, 0)]),
		('D', [(0, 0), (0, .1), (1, .2), (2, .3), (3, .4), (4, .5), (5, .6), (6, .7), (6, 0)]),
		('E', [(0,0), (1, 0.001), (2, 0.2), (3, .1), (4, 0.001), (5,0.001),(6,0)]),
		('F', [(0,0), (1,0.001), (2,0.001), (3,0.001), (4,0.001), (5, 0.001), (6, .1), (6, 0)])
	]
)
def test_get_points(trajectory_table, key, expected):
	series = trajectory_table.loc[key]

	result = polygon.decompose(series)

	assert result == expected



def test_decompose_correct_split_series():
	series = pandas.Series([0.97, 0.0, 0.97, 0.97, 0.97, 0.87, 0.97])
	expected = [0.97, 0.001, 0.97, 0.97, 0.97, 0.87, 0.97]

	result = polygon._decompose_correct_split_series(series)
	assert result.tolist() == pytest.approx(expected)

@pytest.mark.parametrize(
	"timepoint, previous, expected",
	[
		("A", True, "A"),
		("B", True, "A"),
		("C", True, "B"),
		("C", False, "D"),
		("F", False, "G"),
		("G", False, "G")
	]
)
def test_get_neighbor(timepoint, previous, expected):
	index = pandas.Index(list("ABCDEFG"))

	result = polygon.get_neighbor(index, timepoint, previous)

	assert result == expected


@pytest.mark.parametrize(
	"key, expected",
	[
		("A", 1.05), ("B", 0.6), ("C", 4.2), ("D", 2.4), ("E", .3), ("F", 0.05), ("G", 1.1)
	]
)
def test_shoelace(trajectory_table, key, expected):
	# Keep this for now to make sure the area is being calculated correctly.
	series = trajectory_table.loc[key]
	result = areascore.area_of_series(series)

	assert pytest.approx(result, abs = 0.01) == expected



@pytest.mark.parametrize(
	"data, expected",
	[
		([(0, 0), (1, 1), (2, 0), (3, 0), (4, 0), (5, .1), (6, 0)], [[(0, 0), (1, 1), (2, 0)], [(4, 0), (5, .1), (6, 0)]]),
		([(0, 0), (1, 1), (2, 0)], [[(0, 0), (1, 1), (2, 0)]]),
		([], []),
		([(0, 0), (0, .1), (1, .2), (2, .3), (3, .4), (4, .5), (5, .6), (6, .7), (6, 0)],
		[[(0, 0), (0, .1), (1, .2), (2, .3), (3, .4), (4, .5), (5, .6), (6, .7), (6, 0)]])
	]
)
def test_separate(data, expected):
	result = polygon.separate(data)
	assert result == expected

@pytest.mark.parametrize(
	"data, expected",
	[
		([0,1,0,0,.1,0], [[(0,0), (1,1), (2,0.001)], [(3,0.001), (4,.1), (5,0)]]),
		([0,1,0,0,0,.1,0], [[(0,0), (1,1), (2,0.001)], [(4,0.001), (5,0.1), (6,0)]]),
	]
)
def test_separate_again(data, expected):
	series = pandas.Series(data)
	points = polygon.decompose(series)
	result = polygon.separate(points)

	assert result == expected

@pytest.mark.parametrize(
	"data, expected",
	[
		([0,1,0,0,.1,0], [(0,0), (1,1), (2,0.001), (3,0.001), (4,.1), (5,0)]),
		([0,1,0,0,0,.1,0], [(0,0), (1,1), (2,0.001), (3,0.001),(4,0.001), (5,0.1), (6,0)]),
	]
)
def test_decompose(data, expected):
	series = pandas.Series(data)
	result = polygon.decompose(series)
	assert result == expected
