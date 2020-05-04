"""
	Benchmarks the distance calculator. Basically provides an idea of how long it will take to calculate the distance across all pairwise
	trajectories based on how many trajectories there are and how many timepoints were sampled.
"""
import time
from pathlib import Path
from typing import *

import matplotlib.pyplot as plt

plt.style.use("fivethirtyeight")
import pandas
from loguru import logger

from muller import dataio, widgets
from muller.clustering import DistanceCalculator


def load_datasets(filenames: List[Path]) -> List[pandas.DataFrame]:
	tables = list()
	for filename in filenames:
		logger.debug(f"Loading {filename}...")
		t = dataio.import_table(filename, sheet_name = 'trajectory', index = "Trajectory")

		# Remove any unnecessary columns so that the entire table is numeric
		t = t[widgets.get_numeric_columns(t.columns)]

		tables.append(t)
	return tables


def load_highresolution_datasets(filename:Path) -> List[pandas.DataFrame]:

	# Basically want to use this large dataset to generate tables of various sizes.
	table = dataio.import_table(filename, index = "Trajectory")
	table = table[widgets.get_numeric_columns(table.columns)]
	table.columns = list(range(len(table.columns)))
	dfs = list()
	for i in range(2, 10):
		size = int(2 ** i)
		df = table.iloc[:size]
		dfs.append(df)

	return dfs


def get_datasets() -> List[Path]:
	# Use the same datasets as `tests`. May also want to link some larger datasets manually.
	folder_data = Path(__file__).parent.parent.parent / "tests" / "data"

	filenames = [
		folder_data / "tables" / "generic.genotypes.3.xlsx",
		folder_data / "tables" / "generic.genotypes.5.xlsx",
		folder_data / "tables" / "generic.genotypes.10.xlsx",
		# folder_data / "tables" / "real.RMS1-G02.xlsx",
		folder_data / "tables" / "real.nature12344-s2.BYB1-G07.xlsx",
		folder_data / "tables" / "real.B1_muller_try1.xlsx",
		folder_data / "tables_input_trajectories" / "B1_Muller.xlsx",
		folder_data / "tables_input_trajectories" / "P1_Final_Muller.csv",
		folder_data / "tables_input_trajectories" / "P3_Muller.csv"
	]

	return filenames


def benchmark_serial(datasets: List[pandas.DataFrame], threads: int = 1):
	"""
		Benchmarks the serial (single-threaded) distance calculator.
	"""

	calculator = DistanceCalculator(0.0, 0.97, 'binomial', threads = threads)
	benchmarktable = list()
	for dataset in datasets:
		total_elements = len(dataset)
		total_timepoints = len(widgets.get_numeric_columns(dataset.columns))
		start = time.time()
		calculator.run(dataset)
		duration = time.time() - start

		logger.info(f"It took {duration:.2f} seconds to compute {total_elements} elements and {total_timepoints} timepoints.")
		row = {
			'totalTrajectories': total_elements,
			'totalTimepoints':   total_timepoints,
			'duration':          duration,
			'threads':           threads
		}
		benchmarktable.append(row)

	return pandas.DataFrame(benchmarktable)


def plot_banchmarks(data: pandas.DataFrame):
	fig, ax = plt.subplots(figsize = (12, 10))

	ax.scatter("totalTrajectories", "duration", data = data)
	plt.show()


def main():
	filenames_datasets = get_datasets()

	tables_datasets = load_datasets(filenames_datasets)
	#tables_datasets += load_highresolution_datasets()
	dfs = list()
	for i in range(1, 8):
		df_serial = benchmark_serial(tables_datasets, threads = i)
		dfs.append(df_serial)
	d = pandas.concat(dfs)
	d.to_csv("benchmarktable.tsv", sep = "\t")
	plot_banchmarks(dfs[1])


if __name__ == "__main__":
	main()
