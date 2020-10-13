from pathlib import Path

import pandas


class OutputStructure:
	def __init__(self, folder: Path):
		self.filename_genotypes = folder / "tables" / "lineage.trajectory.genotypes.tsv"
		self.filename_edges = folder / "tables" / "lineage.trajectory.lineage.edges.tsv"
		self.filename_figure_mullerpanel = folder / "lineage" / "lineage.trajectory.mullerpanel.png"
		self.filename_muller_table = folder / "tables" / "lineage.trajectory.lineage.muller.tsv"

		self.table_genotypes = pandas.read_csv(self.filename_genotypes, sep = "\t")
		self.table_edges = pandas.read_csv(self.filename_edges, sep = "\t")
		self.table_muller = pandas.read_csv(self.filename_muller_table, sep = "\t")
		# TODO make sure the genotypes table is saved with a Genotype column.