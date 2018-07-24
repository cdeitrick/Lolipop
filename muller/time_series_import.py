from pathlib import Path
import pandas
from typing import List, Tuple

def import_timeseries(filename:Path, timeseries:List = None, sheet_name = 'Sheet1')->Tuple[pandas.DataFrame, pandas.DataFrame]:
	"""
		Reads an excel or csv file. Assumes that the file has the following columns:
	Parameters
	----------
	filename
	timeseries
	sheet_name

	Returns
	-------
	A timeseries dataframe
	the rows correspond to Trajectory, Position, Chromosome and the frequencies for each timepoint.
	"""
	if filename.suffix in {'.xls', '.xlsx'}:
		data = pandas.read_excel(str(filename), sheet_name = sheet_name)
	else:
		data = pandas.read_table(str(filename))

	frequency_columns = list()
	for column in data.columns:
		try:
			int(column)
			frequency_columns.append(column)
		except ValueError:
			pass

	timeseries = data[['Population', 'Trajectory', 'Position'] + frequency_columns]
	#timeseries = timeseries.transpose()

	info = data[['Population', 'Class', 'Mutation' ]]

	return timeseries, info

if __name__ == "__main__":
	filename = "/home/cld100/Documents/github/muller_diagrams/MatLab Muller diagram scripts/B1_muller_try1.xlsx"
	filename = Path(filename)
	import_timeseries(filename, [0,17,25,44,66,75,90])

