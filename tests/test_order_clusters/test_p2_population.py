import pandas
import pytest

from dataio import import_table
from inheritance import order


@pytest.fixture
def p2() -> pandas.DataFrame:
	string = """
		Genotype	0	1	3	4	6	7	9	10	12
		genotype-71	0.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000	1.000
		genotype-65	0.000	0.000	0.000	0.000	0.251	1.000	1.000	1.000	1.000
		genotype-43	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.349	0.000
		genotype-6	0.000	0.000	0.000	0.000	0.175	0.000	0.000	0.000	0.000
		genotype-23	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.162
		genotype-16	0.000	0.058	0.000	0.000	0.000	0.000	0.000	0.000	0.000
		genotype-2	0.000	0.000	0.000	0.000	0.052	0.000	0.000	0.000	0.000
		genotype-21	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.060	0.000
		genotype-39	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.055	0.000
		genotype-8	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.062
		genotype-7	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.054
	"""
	return import_table(string, index = 'Genotype')


@pytest.fixture
def options():
	class Options:
		additive_background_single_cutoff = 1.15
		additive_background_double_cutoff = 1.03
		subtractive_background_single_cutoff = 0.15
		subtractive_background_double_cutoff = 0.03
		derivative_check_cutoff = 0.01
		derivative_detection_cutoff = 0.03
		new_background_base_cutoff = 1.03
		new_background_significant_cutoff = 1.15

	return Options()


def test_p2_population(p2, options):
	expected = {
		'genotype-71': 'genotype-0',
		'genotype-65': 'genotype-71',
		'genotype-43': 'genotype-65',
		'genotype-6':  'genotype-71',
		'genotype-23': 'genotype-65',
		'genotype-16': 'genotype-71',
		'genotype-2':  'genotype-6',
		'genotype-21': 'genotype-43',
		'genotype-39': 'genotype-43',
		'genotype-8':  'genotype-23',
		'genotype-7':  'genotype-8'
	}
	result = order.order_clusters(p2, options)

	assert expected == result.to_dict()
