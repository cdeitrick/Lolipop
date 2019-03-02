import pandas
import pytest

from clustering.metrics import PairwiseCalculationCache


@pytest.fixture
def empty_cache():
	return PairwiseCalculationCache()


@pytest.fixture
def small_cache():
	elements = """
	1	2	.5
	1	3	.6
	1	4	.7
	2	3	.2
	2	4	.3
	3	4	.8
	"""
	string = [i.strip() for i in elements.split("\n")]
	string = "\n".join(i for i in string if i)

	class FakePath:
		def read_text(self) -> str:
			return string

	return PairwiseCalculationCache.read(FakePath())


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
