
import pandas
from dataio import import_table
from clustering.methods.hierarchical_method import hierarchical_method
from clustering.metrics import calculate_pairwise_metric, PairwiseCalculationCache

def bio1() -> pandas.DataFrame:
	string = """
		Trajectory	0	3	4	6	7	9	10	12
		1	0.000	0.653	1.000	1.000	1.000	0.910	0.907	1.000
		2	0.000	0.000	0.646	0.777	0.890	0.512	0.135	0.546
		3	0.000	0.054	0.190	0.131	0.000	0.000	0.053	0.124
		4	0.000	0.080	0.000	0.000	0.000	0.000	0.000	0.226
		5	0.000	0.187	0.000	0.000	0.000	0.000	0.117	0.000
		6	0.000	0.000	0.000	0.000	0.000	0.000	0.759	0.000
	"""
	return import_table(string, index = 'Trajectory')

if __name__ == "__main__":
	expected = [
		["1"],
		["2"],
		["3"],
		["4"],
		["5"],
		["6"]
	]
	print(bio1().to_string())
	print(bio1().info())
	pair_array = calculate_pairwise_metric(bio1(), .03, 0.97, 'binomial')
	result, _ = hierarchical_method(PairwiseCalculationCache(pair_array), similarity_cutoff = 0.03, cluster_method = 'distance')
	from pprint import pprint
	pprint(result)