import pandas
from typing import List, Dict, Tuple

def parse_tree(edges: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Determines the clade and distance from root of every leaf and node in the ggmuller edges table.
	Parameters
	----------
	edges

	Returns
	-------
	pandas.DataFrame
		- index: str
			The leaf/node genotype label
		- columns
			- 'clade': str
				The root genotype that forms the base of the clade.
			- ''distance': int
				The distance from the node/leaf to the root genotype.
	"""
	edges = edges.copy(deep = True)  # To prevent unintended alterations
	leaf_table = edges.set_index('Identity')['Parent']
	clades, iterations = zip(*[determine_clade(leaf_table, i) for i in edges['Identity'].values])
	edges['clade'] = clades
	edges['iterations'] = iterations
	edges = edges.set_index('Identity')
	return edges.sort_values(by = ['clade', 'iterations'])


def determine_clade(parents: pandas.Series, label: str) -> Tuple[str, int]:
	iteration = 0
	while True:
		iteration += 1
		index = parents[label]
		if index != 'genotype-0':
			label = index
		else:
			break
	return label, iteration


def common_elements(left, right):
	common = set(left) & set(right)
	return len(common)


def tokenize(elements: List[str]):
	string = " ".join(elements)
	tokens = string.split(" ")
	return [i.strip() for i in tokens if i]


def group_clades(clade_annotations: Dict[str, List[str]]) -> List[List[str]]:
	"""
		Groups clades by similarity.
		Example
		-------
		clades = {
			'genotype-10': ['dltB '],
			'genotype-11': ['dltB Q111*'],
			'genotype-16': ['PROKKA_00173|PROKKA_00174 '],
			'genotype-6': ['SPAR113_1988 G119C'],
			'genotype-9': ['dltB ']
    	}
    	result = group_clades(clades)
    	result == [['genotype-10', 'genotype-11', 'genotype-9'], ['genotype-16'], ['genotype-6']]
	"""
	tokens = {k: tokenize(v) for k, v in clade_annotations.items()}

	result = list()
	seen = set()
	for key, value in tokens.items():
		if key in seen: continue
		seen.add(key)
		related = [k for k, v in tokens.items() if common_elements(value, v)]
		result.append(related)
		seen = seen | set(related)

	return result


