import pandas
import pytest

from dataio import import_table
from muller.inheritance import order


@pytest.fixture
def three_genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	4	5
	genotype-C	0	0	0	0.3	0.7	1
	genotype-A	0	0	0	0.1	0.5	0.5
	genotype-B	0	0.1	0.15	0.03	0	0
	"""

	return import_table(string, index = 'Genotype')


@pytest.fixture
def five_genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	5	8	13	21	34
	genotype-C	0	0	0	0.3	0.7	1	1	1	1
	genotype-A	0	0	0	0	0.45	0.5	0.55	0.7	0.85
	genotype-E	0	0	0	0	0	0	0.05	0.55	0.5
	genotype-B	0	0.07	0.1	0.02	0.01	0	0	0	0
	genotype-D	0	0	0	0	0	0	0.07	0	0.01
	"""
	return import_table(string, index = 'Genotype')


@pytest.fixture
def options():
	class Options:
		additive_background_single_cutoff = 0.15
		additive_background_double_cutoff = 0.03
		subtractive_background_single_cutoff = 0.15
		subtractive_background_double_cutoff = 0.03
		derivative_check_cutoff = 0.01
		derivative_detection_cutoff = 0.03
		new_background_base_cutoff = 1.03
		new_background_significant_cutoff = 1.15

	return Options()


def test_three_genotypes(three_genotypes, options):
	expected = {
		'genotype-A': 'genotype-C',
		'genotype-B': 'genotype-0',
		'genotype-C': 'genotype-0'
	}
	clusters = order.order_clusters(three_genotypes, options)

	assert expected == clusters.to_dict()


def test_five_genotypes(five_genotypes, options):
	expected = {
		'genotype-A': 'genotype-C',
		'genotype-B': 'genotype-0',
		'genotype-C': 'genotype-0',
		'genotype-D': 'genotype-A',
		'genotype-E': 'genotype-A'
	}
	result = order.order_clusters(five_genotypes, options)
	assert expected == result.to_dict()
