from pathlib import Path
import pandas
from typing import Tuple


def import_timeseries(filename: Path, sheet_name = 'Sheet1') -> Tuple[pandas.DataFrame, pandas.DataFrame]:
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
	if filename.suffix in {'.xls', '.xlsx'}:
		data = pandas.read_excel(str(filename), sheet_name = sheet_name)
	else:
		data = pandas.read_table(str(filename))

	# Extract the columns which indicate timepoints of observations. Should be integers.
	frequency_columns = list()
	for column in data.columns:
		try:
			int(column)
			frequency_columns.append(column)
		except ValueError:
			pass

	# Extract the columns with the trajectory identifiers and frequencies at each timepoint.
	timeseries = data[['Population', 'Trajectory', 'Position'] + frequency_columns]
	# timeseries = timeseries.transpose()

	# Extract metadata for each trajectory.
	info = data[['Population', 'Class', 'Mutation']]

	return timeseries, info


if __name__ == "__main__":
	pass