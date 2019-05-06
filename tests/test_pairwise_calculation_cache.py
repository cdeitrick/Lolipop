import pandas
import pytest

from clustering.metrics import PairwiseCalculationCache


@pytest.fixture
def empty_cache():
	return PairwiseCalculationCache()


@pytest.fixture
def small_cache():
	data = {
		('1', '2'): .5,
		('1', '3'): .6,
		('1', '4'): .7,
		('2', '3'): .2,
		('2', '4'): .3,
		('3', '4'): .8
	}

	return PairwiseCalculationCache(data)


def test_cache_empty(empty_cache, small_cache):
	assert bool(empty_cache) == False
	assert bool(small_cache) == True  # For clarity.


def test_read(small_cache):
	expected = {
		('1', '2'): .5,
		('2', '1'): .5,

		('1', '3'): .6,
		('3', '1'): .6,

		('1', '4'): .7,
		('4', '1'): .7,

		('2', '3'): .2,
		('3', '2'): .2,

		('2', '4'): .3,
		('4', '2'): .3,

		('3', '4'): .8,
		('4', '3'): .8
	}

	assert expected == small_cache.pairwise_values


def test_squareform(small_cache):
	expected = [
		[0, .5, .6, .7],
		[.5, 0, .2, .3],
		[.6, .2, 0, .8],
		[.7, .3, .8, 0]
	]
	labels = "1 2 3 4".split()
	expected_df = pandas.DataFrame(expected, columns = labels, index = labels)
	pandas.testing.assert_frame_equal(expected_df, small_cache.squareform())


def test_asdict(small_cache):
	expected = {
		('1', '2'): .5,
		('1', '3'): .6,
		('1', '4'): .7,
		('2', '1'): .5,
		('2', '3'): .2,
		('2', '4'): .3,
		('3', '1'): .6,
		('3', '2'): .2,
		('3', '4'): .8,
		('4', '1'): .7,
		('4', '2'): .3,
		('4', '3'): .8
	}

	assert small_cache.asdict() == expected


def test_unique(small_cache):
	expected = [
		('1', '2'), ('1', '3'), ('1', '4'),
		('2', '3'), ('2', '4'),
		('3', '4')
	]

	result = sorted(small_cache.unique())
	assert result == expected


def test_update(small_cache):
	new_elements = {
		('15', '16'): 1,
		('14', '1'):  2
	}

	small_cache.update(new_elements)

	assert small_cache.get('14', '1') == 2
	assert small_cache.get('1', '14') == 2
