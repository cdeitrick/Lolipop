from pathlib import Path
from typing import *
import pandas
import matplotlib.pyplot as plt
from muller import treetools
def generate_lineage_table(genotypes:pandas.DataFrame, lineages:Dict[str,str]):


	groups = genotypes.groupby(by = lambda s: lineages[s])
	lineage_table = list()
	for lineage, group in groups:
		average = group.median()
		lineage_table.append(average)
	return pandas.DataFrame(lineage_table)



def check_lineage(genotypes:pandas.DataFrame, lineages:pandas.Series):

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


def main():
	filename = "/home/cld100/Documents/github/muller_diagrams/manuscript_scripts/03 - similarity cutoff comparison/results/sc-1/real.nature12344-s2.BYB1-G07.genotypes.tsv"
	filename_edges = "/home/cld100/Documents/github/muller_diagrams/manuscript_scripts/03 - similarity cutoff comparison/results/sc-1/tables/real.nature12344-s2.BYB1-G07.edges.tsv"
	table = pandas.read_csv(filename, sep = "\t")
	table = table.set_index('Genotype')
	if 'members' in table.columns:
		table.pop('members')
	edges = pandas.read_csv(filename_edges, sep = "\t")
	edges = edges.set_index('Identity')
	lineages = dict()
	for identity in edges.index:
		lineage = get_lineage(edges, identity)
		lineages[identity] = lineage

	print(table.to_string())
	lineage_table = generate_lineage_table(table, lineages)
	table['lineage'] = [get_lineage(edges, i) for i in table.index]
	print(table.to_string())
	print(lineage_table.to_string())
	print(lineage_table.sum()-1)


if __name__ == "__main__":
	main()