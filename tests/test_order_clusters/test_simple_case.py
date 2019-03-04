import pandas
import pytest

from dataio import import_table


@pytest.fixture
def three_genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	4	5
	Trajectory-C	0	0	0	0.3	0.7	1
	Trajectory-A	0	0	0	0.1	0.5	0.5
	Trajectory-B	0	0.1	0.15	0.03	0	0
	"""

	return import_table(string, index = 'Genotype')


@pytest.fixture
def options():
	class Options:
		additive_check_single_cutoff = 0.15
		additive_check_double_cutoff = 0.03
		subtractive_check_single_cutff = 0.15
		subtractive_check_double_cutoff = 0.03
		derivative_check_cutoff = 0.01
		derivative_detection_cutoff = 0.03

	return Options()


def test_three_genotypes(three_genotypes):
	pass
