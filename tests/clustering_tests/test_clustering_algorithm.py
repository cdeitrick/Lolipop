from pathlib import Path
import pytest
import pandas
from loguru import logger
from muller import clustering
from typing import List, Dict
@pytest.fixture
def trajectories()->pandas.DataFrame:
	folder = Path.home() / "storage" / "projects" / "muller"
	filename_input = folder / "High-resolution lineage" /"source data figure 2"/ "source data figure 2.edited.YPD.tsv"

	table = pandas.read_csv(filename_input, sep = "\t").set_index('Trajectory')
	logger.debug(list(table.columns))
	table = table.drop('Population', axis = 'columns')
	print(list(table.columns))

	return table

def _write_genotypes_to_file(filename:Path, genotypes:Dict[str,List[str]]):
	with filename.open('w') as file1:
		for k,v in genotypes.items():
			line = "|".join(v)
			file1.write(f"{k}\t{line}\n")

def test_clustering_algorithm(trajectories):


	# (self, metric: str, dlimit: float, slimit:float, flimit: float, pvalue: float,
	# 			starting_genotypes: Optional[List[List[str]]] = None, threads: Optional[int] = None):

	workflow = clustering.ClusterMutations(
		metric = "binomial",
		dlimit = 0.01,
		slimit = 0.15,
		flimit = 0.99,
		pvalue = 0.05,
	)

	result = workflow.run(trajectories)

	f = result.genotype_members
	folder = Path.home() / "storage" / "projects" / "muller"
	filename_output = folder / "High-resolution lineage" / "source data figure 2" / "source data figure 2.edited.YPD.genotypes.tsv"


	_write_genotypes_to_file(filename_output, f)
	print(f)