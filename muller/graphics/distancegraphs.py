from pathlib import Path
from typing import List, Optional

from loguru import logger

from muller import clustering, dataio
import seaborn
import matplotlib.pyplot as plt
plt.style.use("/home/cld100/Documents/sandbox/matplotlibrc")
def generate_distance_plot(distances: List[float], similarity_cutoff: float, filename:Optional[Path] = None, ax:plt.Axes = None):
	""" Shows the spread of computed pairwise distances."""
	if not ax:
		fig, ax = plt.subplots(figsize = (12, 10))

	seaborn.distplot(distances, ax = ax, kde = False, rug = True, bins = 20)
	ax.axvline(similarity_cutoff, color = 'red')

	ax.set_title("Pairwise distances between each pair of trajectories")
	ax.set_xlabel("Distance")
	ax.set_ylabel("Count")
	ax.set_xlim(0, max(distances))
	plt.tight_layout()
	if filename:
		plt.savefig(filename)
	else:
		plt.show()


def main():
	folder_github = Path.home() / "Documents" / "github" / "muller_diagrams"
	folder_project = Path.home() / "storage" / "projects" / "muller"
	folder_data = folder_github / "tests" / "data"
	filename_b1 = folder_data / "tables" / "real.B1_muller_try1.edited.xlsx"
	filename_5 = folder_data / "tables" / "generic.genotypes.5.xlsx"
	filename_high = folder_project / "High-resolution lineage" / "source data figure 2" / "source data figure 2.edited.YPD.tsv"

	table = dataio.import_table(filename_b1, index = 'Trajectory', sheet_name = 'trajectory')

	table = table[[i for i in table.columns if isinstance(i, int)]]
	workflow = clustering.ClusterMutations(
		metric = 'binomial',
		dlimit = 0.03,
		pvalue = 0.05,
		slimit = 0.15,
		flimit = 0.97
	)
	workflow_p = clustering.ClusterMutations(
		metric = 'binomialp',
		dlimit = 0.03,
		pvalue = 0.05,
		slimit = 0.15,
		flimit = 0.97
	)

	result = workflow.run(table, distance_cutoff = 0.05)
	resultp = workflow_p.run(table, distance_cutoff = 0.05)

	#values = list(result.matrix_distance.pairwise_values.values())

	result_left = list(result.matrix_distance.pairwise_values.values())
	result_right = list(resultp.matrix_distance.pairwise_values.values())

	ax = seaborn.scatterplot(result_left, result_right)
	ax.axvline(result.clusterdata.similarity_cutoff, color = 'red')
	ax.axhline(resultp.clusterdata.similarity_cutoff, color = 'green')
	ax.set_xlabel("binomial")
	ax.set_ylabel("binomialp")

	plt.show()

	from pprint import pprint
	# for k,v in sorted(result.matrix_distance.pairwise_values.items()):
	#	print(k,v)
	#pprint(result.genotype_members)

	#generate_distance_plot(values, similarity_cutoff = result.clusterdata.similarity_cutoff)


if __name__ == "__main__":
	main()
