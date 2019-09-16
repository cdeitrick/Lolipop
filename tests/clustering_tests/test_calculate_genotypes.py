from io import StringIO

import pandas
import pytest

from muller.clustering import ClusterMutations
from muller.dataio import import_table

trajectory_csv = "Trajectory,0,17,25,44,66,75,90\n" \
				 "1,0,0.0,0.261,1.0,1.0,1.0,1.0\n" \
				 "2,0,0.0,0.0,0.525,0.454,0.911,0.91\n" \
				 "3,0,0.0,0.0,0.147,0.45,0.924,0.887\n" \
				 "4,0,0.0,0.0,0.0,0.211,0.811,0.813\n" \
				 "6,0,0.0,0.0,0.0,0.0,1.0,1.0\n" \
				 "7,0,0.0,0.0,0.273,0.781,1.0,1.0\n" \
				 "8,0,0.0,0.0,0.0,0.345,0.833,0.793\n" \
				 "10,0,0.0,0.117,0.0,0.0,0.0,0.103\n" \
				 "11,0,0.0,0.0,0.108,0.151,0.0,0.0\n" \
				 "13,0,0.0,0.0,0.0,0.258,0.057,0.075\n" \
				 "14,0,0.38,0.432,0.0,0.0,0.0,0.0\n" \
				 "16,0,0.0,0.0,0.0,0.209,0.209,0.0\n" \
				 "20,0,0.0,0.0,0.138,0.295,0.0,0.081\n"


@pytest.fixture
def genotype_generator() -> ClusterMutations:
	generator = ClusterMutations(
		method = 'hierarchy',
		metric = 'binomial',
		dlimit = 0.03,
		flimit = 0.97,
		pvalue = 0.05,
		dbreakpoint = 0.10,
		breakpoints = [],
		starting_genotypes = [],
		trajectory_filter = None
	)
	return generator


def test_calculate_mean_genotype(genotype_generator):
	test_genotypes = [
		['7'], ['4', '8'], ['3', '2'], ['13', '20', '11']
	]
	trajectories = pandas.read_csv(StringIO(trajectory_csv))
	trajectories['Trajectory'] = trajectories['Trajectory'].astype(str)
	trajectories = trajectories.set_index('Trajectory')

	expected_csv = """
		Genotype	0	17	25	44	66	75	90	members
		genotype-1	0.0	0.0	0.0	0.273	0.781	1.0	1.0	7
		genotype-2	0	0	0	0	0.278	0.822	0.803	4|8
		genotype-3	0	0	0	0.336	0.452	0.9175	0.8985	3|2
		genotype-4	0	0	0	0.082	0.234666666666667	0.019	0.052	13|20|11
		"""
	expected_mean = import_table(expected_csv, index = 'Genotype')
	output = genotype_generator.calculate_mean_genotype(test_genotypes, trajectories)

	pandas.testing.assert_frame_equal(expected_mean, output)
