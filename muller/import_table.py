from pathlib import Path
import pandas
from typing import Tuple, List


def get_numeric_columns(columns: List[str]) -> List[str]:
	numeric_columns = list()
	for column in columns:
		try:
			int(column)
			numeric_columns.append(column)
		except ValueError:
			pass
	return numeric_columns


def import_table(filename: Path, sheet_name:str) -> pandas.DataFrame:
	if filename.suffix in {'.xls', '.xlsx'}:

		data: pandas.DataFrame = pandas.read_excel(str(filename), sheet_name = sheet_name)

	else:
		sep = '\t' if filename.suffix in {'.tsv', '.tab'} else ','
		data: pandas.DataFrame = pandas.read_table(str(filename), sep = sep)
	return data


def import_genotype_table(filename: Path) -> pandas.DataFrame:
	data = import_table(filename)

	if 'Genotype' in data.columns:
		key_column = 'Genotype'
	elif 'Unnamed: 0' in data.columns:
		key_column = 'Unnamed: 0'
	else:
		message = "One of the columns needs to be labeled 'Genotype'"
		raise ValueError(message)

	if 'members' not in data.columns:
		message = "The genotype must have a 'members' column with the names of all trajectories contained in the genotype. Individual trajectory names must be separated by '|'"
		raise ValueError(message)

	data = data[[key_column, 'members'] + get_numeric_columns(data.columns)]

	data = data.set_index(key_column)
	data.index.name = 'Genotype'

	def _convert_to_numeric(col):
		try:
			return int(col)
		except:
			return col

	data = data.rename(_convert_to_numeric, axis = 'columns')
	return data


def import_trajectory_table(filename: Path, sheet_name = 'Sheet1') -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	"""
		Reads an excel or csv file. Assumes that the file has the following columns:
		- Population:
		- Trajectory:
		- Position:
	Parameters
	----------
	filename: Path
		The table containing the trajectories and associated metadata. Can be an excel sheet or comma/tab delimited file.
	sheet_name: str; Default 'Sheet1'
		Indicates which sheet contains the data, if an excel table is given.
	Returns
	-------
	A timeseries dataframe
		- Columns
			- Population: str
				Name of the population. ex. 'B2'
			- Trajectory: int
				Identifies a unique mutation based on population and posiiton. Should be sorted starting from 1
			- Position: int
				Position of the mutation.
			* timeseries
				The timeseries points will correspond to the timepoints included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	A dataframe with metadata for each trajectory. Includes Population, Class (ex 'SNP'), and Mutation ('C>T')
	"""

	# Read in the data table.
	data = import_table(filename, sheet_name)

	key_column = 'Trajectory'
	# Extract the columns which indicate timepoints of observations. Should be integers.
	frequency_columns = get_numeric_columns(data.columns)

	# Extract the columns with the trajectory identifiers and frequencies at each timepoint.

	data = data.set_index(key_column)
	timeseries = data[frequency_columns]
	for column in frequency_columns:
		if max(timeseries[column]) > 1.0:
			timeseries[column] /= 100


	# Extract metadata for each trajectory.
	try:
		potential_columns = ['Population','Position', 'Class', 'Gene', 'Mutation', 'Annotation']
		info = data[[i for i in potential_columns if i in data.columns]]
		#info = data[potential_columns]
	except:
		info = None
	timeseries = timeseries.rename(columns = int)
	return timeseries, info


if __name__ == "__main__":
	pass
