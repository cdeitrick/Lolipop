import itertools
from pathlib import Path
from typing import Dict, Tuple

import pandas
import pytest

from muller import dataio
from muller.inheritance import areascore, scoring
from .filenames import FILENAME_TRUTHSET


#################################################################################
############################### Fixtures ########################################
#################################################################################

@pytest.fixture
def full_overlap() -> Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, 0, 0, 0, 1, 1])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


@pytest.fixture
def no_overlap() -> Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, .3, .4, 0, 0, 0])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


@pytest.fixture
def partial_overlap() -> Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, 0, .1, .2, 1, 1])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


@pytest.fixture
def truthset_common_area() -> pandas.DataFrame:
	t = dataio.import_table(FILENAME_TRUTHSET, sheet_name = 'commonarea', index = 'index')
	return t


@pytest.fixture
def truthset_overlap() -> pandas.DataFrame:
	t = dataio.import_table(FILENAME_TRUTHSET, sheet_name = 'overlap', index = 'index')
	t.loc['F'] = t.loc['F'].astype(str)
	return t


@pytest.fixture
def truthset_or_area() -> pandas.DataFrame:
	return dataio.import_table(FILENAME_TRUTHSET, sheet_name = 'orarea', index = 'index')


@pytest.fixture
def truthset_xor_area() -> pandas.DataFrame:
	return dataio.import_table(FILENAME_TRUTHSET, sheet_name = 'xorarea', index = 'index')


@pytest.fixture
def datatables() -> pandas.DataFrame:
	df = pandas.read_excel(FILENAME_TRUTHSET, sheet_name = 'data')
	df = df.set_index('Trajectory')
	return df


@pytest.fixture
def actualarea() -> Dict[str, float]:
	A = [0.1 / 2, 0.1, 0.1 + (0.1 / 2), 0.2, 0.2 + (0.1 / 2), 0.3]
	B = [0.1] * 6
	C = (0.4 * 6) + (0.6 * 6 / 2)
	D = (.1 * 6) + (.6 * 6 / 2)
	E = 0.1 + (0.2 * 2 / 2)
	F = 0.05
	G = 1 + (2 * 0.2 / 2)
	d = {
		'A': sum(A),
		'B': sum(B),
		'C': C,
		'D': D,
		'E': E,
		'F': F,
		'G': G
	}

	return d


@pytest.fixture
def truthset_subsetscore() -> pandas.DataFrame:
	s = """
		index	A	B	C	D	E	F	G
		A		-2	2	2	0	-2	0
		B	1		2	2	0	-2	0
		C	-2	-2		0			
		D	-2	-2	2		-2	-2	-2
		E	1	1	2	2		0	0
		F	2	2	2	2	0		0
		G	0	0	2	2	0	0	
	"""
	return dataio.import_table(FILENAME_TRUTHSET, sheet_name = 'score.area.expected', index = 'unnested/nested')


@pytest.fixture
def truthset_issubsetbool() -> pandas.DataFrame:
	filename = FILENAME_TRUTHSET
	table = dataio.import_table(filename, index = 'unnested/nested', sheet_name = 'issubsetbool')
	table = table.fillna(True)  # Each series should be a subset of itself.
	table = table.astype(bool)
	return table


#################################################################################
############################### Small Tests #####################################
#################################################################################

def test_calculate_overlapping_area(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	right_area = areascore.area_of_series(right)
	overlap_area = areascore.calculate_common_area(left, right)
	assert overlap_area == right_area

	left, right = no_overlap
	overlap_area = areascore.calculate_common_area(left, right)
	assert overlap_area == 0

	left, right = partial_overlap
	overlap_area = areascore.calculate_common_area(left, right)
	assert overlap_area == right_area

	left = pandas.Series([0.01, 0.279, 0.341, 0.568, 0.708, 0.913, 0.756, 0.455, 0.399, 0.13, 0.041])
	right = pandas.Series([0, 0, 0, 0, 0, 0.247, 0.388, 0.215, 0.399, 0.13, 0.028])
	overlap_area = areascore.calculate_common_area(left, right)
	assert pytest.approx(overlap_area) == areascore.area_of_series(right)


def test_is_subset(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	assert areascore.is_subset(left, right)
	assert False == areascore.is_subset(right, left)

	left, right = partial_overlap
	assert True == areascore.is_subset(left, right)
	assert False == areascore.is_subset(right, left)
	left, right = no_overlap
	assert False == areascore.is_subset(left, right)
	assert False == areascore.is_subset(right, left)


@pytest.mark.parametrize(
	"series, expected",
	[
		([0, 0.2, 0.3, 0.4, 0.5], 1.4),
		([0.3, 0.3, 0.3, 0.3, 0.3], 1.5),
		([0, .1, .1, .2, .2, .3, .3], 1.2),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], 1.029),
		([0, 0, 0, 0, 0], 0),
		([], 0)

	]
)
def test_area_of_series(series, expected):
	s = pandas.Series(series)
	result = areascore.area_of_series(s)

	assert pytest.approx(result, expected)


#################################################################################
############################### Permutation Tests ###############################
#################################################################################

def helper_for_datatables(key, dtable, ttable):
	""" Removes some boilerplate code from the tests."""
	l, r = key
	ls = dtable.loc[l]
	rs = dtable.loc[r]

	e = ttable.at[r, l]

	return ls, rs, e


@pytest.mark.parametrize("key", list(itertools.permutations('ABCDEFG', 2)))  # Use permutations to test if result changes if the parameters switch.
def test_calculate_common_area(datatables, truthset_common_area, key):
	left, right, expected = helper_for_datatables(key, datatables, truthset_common_area)

	result = areascore.calculate_common_area(left, right)

	assert pytest.approx(result, abs = 0.05) == expected


@pytest.mark.parametrize("key", list(itertools.permutations('ABCDEFG', 2)))  # Use permutations to test if result changes if the parameters switch.
def test_X_or_Y(datatables, truthset_or_area, key):
	left, right, expected = helper_for_datatables(key, datatables, truthset_or_area)

	result = areascore.X_or_Y(left, right)
	assert pytest.approx(result, abs = 0.01) == expected


@pytest.mark.parametrize("key", list(itertools.permutations('ABCDEFG', 2)))  # Use permutations to test if result changes if the parameters switch.
def test_X_xor_Y(datatables, truthset_xor_area, key):
	left, right, expected = helper_for_datatables(key, datatables, truthset_xor_area)

	result = areascore.X_xor_Y(left, right)

	assert pytest.approx(result, abs = 0.01) == expected

@pytest.mark.parametrize("key", list(itertools.permutations('ABCDEFG', 2)))  # Use permutations to test if result changes if the parameters switch.
def test_area_score(datatables, truthset_subsetscore, key):
	left, right, expected = helper_for_datatables(key, datatables, truthset_subsetscore)
	scorer = scoring.Score(0.03, 0.97, 0.05)
	result = scorer.calculate_score_area(left, right)

	assert result == expected


@pytest.mark.parametrize("key", list(itertools.permutations('ABCDEFG', 2)))  # Use permutations to test if result changes if the parameters switch.
def test_issubsetbool(datatables, truthset_issubsetbool, key):
	left, right, expected = helper_for_datatables(key, datatables, truthset_issubsetbool)

	result = areascore.is_subset(left, right)
	assert result == expected


def test_problematic_dataset():
	# This dataset is causing a TopologyError when calculating the common area
	string = """
	Genotype	0	1	2	3	4	5
	genotype-C	0	0	0	0.3	0.7	1
	genotype-A	0	0	0	0.1	0.5	0.5
	genotype-B	0	0.1	0.15	0.03	0	0
	"""

	#table = import_table(string, index = 'Genotype')