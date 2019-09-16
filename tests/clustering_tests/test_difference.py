import pandas
import pytest

from muller.clustering.methods import twostep_unlink


@pytest.fixture
def trajectories() -> pandas.Series:
	data = {
		('1', '2'): .5,
		('1', '3'): .7,
		('1', '4'): .2,
		('2', '1'): .5,
		('2', '3'): .7,
		('2', '4'): .1,
		('3', '1'): .7,
		('3', '2'): .7,
		('3', '4'): .01,
		('4', '1'): .2,
		('4', '2'): .1,
		('4', '3'): .01
	}

	return pandas.Series(data)


def test_find_weakest_pair(trajectories):
	assert twostep_unlink._find_weakest_pair(trajectories, .1) == ('2', '4')


def test_divide_genotype(trajectories):
	genotype = ['1', '2', '3', '4']

	expected_genotype1 = ['2', '1', '3']
	expected_genotype2 = ['4']

	result1, result2 = twostep_unlink._divide_genotype(genotype, trajectories, 0.1)

	assert expected_genotype1 == result1
	assert expected_genotype2 == result2


def test_unlink_unrelated_trajectories(trajectories):
	genotypes = [['1', '2', '3', '4']]

	expected = [['4'], ['2', '1', '3']]

	result = twostep_unlink.unlink_unrelated_trajectories(genotypes, trajectories.to_dict(), 0.1)

	assert result == expected
