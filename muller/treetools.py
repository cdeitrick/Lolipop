from typing import Dict, List, Tuple

import pandas


def get_child_nodes(tree: pandas.DataFrame, label: str) -> List[str]:
	"""Retrieves all child node for the given label from the tree. Includes `label` in the output"""

	children = list()
	for index in tree.index:
		parents = get_parent_nodes(tree, index)
		if label in parents:
			children.append(index)
	return children


def get_parent_nodes(tree: pandas.DataFrame, label: str) -> List[str]:
	"""Retrieves all parent nodes for the given node."""
	parents = []
	while label != 'genotype-0':
		try:
			parent = tree.loc[label]['Parent']
		except KeyError:
			parent = 'genotype-0'
		parents.append(parent)
		label = parent
	return parents


def parse_tree(edges: pandas.Series) -> pandas.DataFrame:
	"""
		Determines the clade and distance from root of every leaf and node in the ggmuller edges table.
	Parameters
	----------
	edges

	Returns
	-------
	pandas.Series
		- index: str
			The leaf/node genotype label
		- columns
			- 'clade': str
				The root genotype that forms the base of the clade.
			- ''distance': int
				The distance from the node/leaf to the root genotype.
	"""
	edges = edges.copy(deep = True)  # To prevent unintended alterations

	#lineage_table = edges.set_index('Identity')['Parent']
	lineage_table = edges

	clades, iterations = zip(*[determine_clade(lineage_table, i) for i in lineage_table.index])
	leaf_table = lineage_table.to_frame().reset_index()
	leaf_table['clade'] = clades
	leaf_table['iterations'] = iterations
	leaf_table = leaf_table.set_index('Identity')
	return leaf_table.sort_values(by = ['clade', 'iterations'])


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


def common_elements(left, right) -> int:
	""" Returns the number of elements that are common between `left` and `right`"""
	common = set(left) & set(right)
	return len(common)


def tokenize(elements: List[str]):
	""" Converts a list of strings into a list of unique words found in all the string elements."""
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
