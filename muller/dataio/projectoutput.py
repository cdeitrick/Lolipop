from pathlib import Path
from typing import *
import pandas
class OutputStructure:
	def __init__(self, folder: Path):
		self.filename_genotypes = folder / "tables" / "lineage.genotypes.tsv"
		self.filename_edges = folder / "tables" / "lineage.lineage.edges.tsv"

		self.table_genotypes = pandas.read_csv(self.filename_genotypes, sep = "\t")
		self.table_edges = pandas.read_csv(self.filename_edges, sep = "\t")
