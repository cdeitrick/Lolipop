from typing import *

import matplotlib.pyplot as plt
import pandas


def generate_lineage_table(genotypes:pandas.DataFrame, lineages:Dict[str,str]):


	groups = genotypes.groupby(by = lambda s: lineages[s])
	lineage_table = list()
	for lineage, group in groups:
		average = group.max()
		lineage_table.append(average)
	return pandas.DataFrame(lineage_table)


def calculate_score(table:pandas.DataFrame, flimit:float = 0.97)->float:
	ancestor = get_ancestor_series(table, 0.97)
	genotypes = table[[i for i in table.columns if i not in ancestor.index]]
	frequencies = genotypes.sum()

	quad_frequencies = frequencies**2
	score = quad_frequencies.sum()/len(table.columns)
	print(frequencies)
	print(quad_frequencies)
	print(score)


def check_lineage(genotypes:pandas.DataFrame, lineages:pandas.Series):

	ancestor = get_ancestor_series(genotypes, 0.97)
	genotypes = genotypes[[i for i in genotypes.columns if i not in ancestor.index]]
	frequencies = genotypes.sum()


	fig, ax = plt.subplots(figsize = (10,10))

	ax.plot(frequencies)

	plt.show()


def get_lineage(edges:pandas.Series, identity:str):
	while True:
		parent = edges.loc[identity].iloc[0]
		if parent != 'genotype-0':
			identity = parent
		else:
			break
	return identity

def get_first_nonancestor_timepoint(table:pandas.DataFrame, flimit:float = 0.97)->int:
	""" Returns the first timepoint that the ancestor was not present in."""
	for column_label in table.columns:

		column = table[column_label]
		if column.sum() >= flimit:
			return column_label
	else:
		# If the ancestor was never eradicated, return the last timepoint.
		return column_label

def get_ancestor_series(table:pandas.DataFrame, flimit:float = 0.97)->pandas.Series:

	endpoint = get_first_nonancestor_timepoint(table, flimit)
	early_table = table[[i for i in table.columns if i < endpoint]]
	nonancestor = early_table.sum()
	ancestor = 1 - nonancestor

	return ancestor


def read_tables()->Tuple[pandas.DataFrame, pandas.DataFrame]:
	filename = "/home/cld100/Documents/github/muller_diagrams/manuscript_scripts/03 - similarity cutoff comparison/results/sc-1/real.nature12344-s2.BYB1-G07.genotypes.tsv"
	filename_edges = "/home/cld100/Documents/github/muller_diagrams/manuscript_scripts/03 - similarity cutoff comparison/results/sc-1/tables/real.nature12344-s2.BYB1-G07.edges.tsv"
	table = pandas.read_csv(filename, sep = "\t")
	table = table.set_index('Genotype')
	if 'members' in table.columns:
		table.pop('members')
	edges = pandas.read_csv(filename_edges, sep = "\t")
	edges = edges.set_index('Identity')

	table.columns = [int(i) for i in table.columns]

	return table, edges

def main():
	table_genotypes, edges = read_tables()
	lineages = dict()
	for identity in edges.index:
		lineage = get_lineage(edges, identity)
		lineages[identity] = lineage

	lineage_table = generate_lineage_table(table_genotypes, lineages)
	#table_genotypes['lineage'] = [get_lineage(edges, i) for i in table_genotypes.index]
	calculate_score(lineage_table)

if __name__ == "__main__":
	main()