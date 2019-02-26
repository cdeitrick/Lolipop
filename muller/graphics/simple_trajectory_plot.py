"""
	Designed to plot simple tables of trajectories which are named based on expected genotype.
	Ex.
	Trajectory	0	1	2	3	5	8	13	21	34
	trajectory-A1	0	0	0	0.1	0.5	0.5	0.61	0.7	0.8
	trajectory-A2	0	0	0	0.06	0.35	0.4	0.65	0.75	0.8
	trajectory-A3	0	0	0	0	0.45	0.5	0.55	0.7	0.85
	trajectory-B1	0	0.1	0.15	0.03	0	0	0	0	0
	trajectory-B2	0	0.07	0.1	0.02	0.01	0	0	0	0
	trajectory-C1	0	0	0	0.3	0.7	1	1	1	1
	trajectory-D1	0	0	0	0	0	0	0.1	0	0
	trajectory-D2	0	0	0	0	0	0	0.07	0	0.01
	trajectory-E1	0	0	0	0	0	0	0.1	0.5	0.5
	trajectory-E2	0	0	0	0	0	0	0.05	0.55	0.5
	trajectory-E3	0	0	0	0	0	0	0	0.6	0.5
"""
from pathlib import Path
import pandas
from typing import List, Dict
import matplotlib.pyplot as plt
import palette
def sort_into_genotypes(labels:List[str])->Dict[str,List[str]]:
	data = dict()
	for label in labels:
		trajectory_id = label.split('-')[-1]
		genotype_id = trajectory_id[0]
		if genotype_id in data:
			data[genotype_id].append(label)
		else:
			data[genotype_id] = [label]
	return data
def simple_plot(table:pandas.DataFrame, output_filename:Path):

	genotypes = sort_into_genotypes(table.index)

	fig, ax = plt.subplots(figsize = (20,10))
	for (genotype_id, member_trajectories), color in zip(genotypes.items(), palette.DISTINCTIVE_PALETTE):
		group = table.loc[member_trajectories]

		for index, row in group.iterrows():
			ax.plot(row.index, row.values, color = color)
	plt.savefig(str(output_filename))

if __name__ == "__main__":
	from import_data import import_table_from_string
	string = """
	Trajectory	0	1	2	3	5	8	13	21	34
	trajectory-A1	0	0	0	0.1	0.5	0.5	0.61	0.7	0.8
	trajectory-A2	0	0	0	0.06	0.35	0.4	0.65	0.75	0.8
	trajectory-A3	0	0	0	0	0.45	0.5	0.55	0.7	0.85
	trajectory-B1	0	0.1	0.15	0.03	0	0	0	0	0
	trajectory-B2	0	0.07	0.1	0.02	0.01	0	0	0	0
	trajectory-C1	0	0	0	0.3	0.7	1	1	1	1
	trajectory-D1	0	0	0	0	0	0	0.1	0	0
	trajectory-D2	0	0	0	0	0	0	0.07	0	0.01
	trajectory-E1	0	0	0	0	0	0	0.1	0.5	0.5
	trajectory-E2	0	0	0	0	0	0	0.05	0.55	0.5
	trajectory-E3	0	0	0	0	0	0	0	0.6	0.5
	"""
	table = import_table_from_string(string, index = 'Trajectory')
	simple_plot(table, "five_genotypes.png")

