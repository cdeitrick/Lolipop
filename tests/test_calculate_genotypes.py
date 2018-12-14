import unittest
from io import StringIO

from muller_genotypes.calculate_genotypes import *
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


class TestCalculateGenotypes(unittest.TestCase):
	def test_calculate_mean_genotype(self):
		test_genotypes = [
			['7'], ['4', '8'], ['3', '2'], ['13', '20', '11']
		]
		trajectories = pandas.read_csv(StringIO(trajectory_csv))
		trajectories['Trajectory'] = trajectories['Trajectory'].astype(str)

		trajectories.set_index('Trajectory', inplace = True)

		expected_csv = "Genotype	0	17	25	44	66	75	90	members\n" \
					   "genotype-1	0	0.0	0.0	0.273	0.781	1.0	1.0	7\n" \
					   "genotype-2	0	0	0	0	0.278	0.822	0.803	4|8\n" \
					   "genotype-3	0	0	0	0.336	0.452	0.9175	0.8985	3|2\n" \
					   "genotype-4	0	0	0	0.082	0.234666666666667	0.019	0.052	13|20|11"
		expected_mean = pandas.read_csv(StringIO(expected_csv), sep = '\t', index_col = 'Unnamed:0')
		expected_mean['0'] = expected_mean['0'].astype(float)

		output = calculate_mean_genotype(test_genotypes, trajectories)

		pandas.testing.assert_frame_equal(expected_mean, output)


