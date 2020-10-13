from pathlib import Path
from typing import *

import pandas


def read_edges_filename(path:Path)->pandas.DataFrame:
	# Use a dataframe in case the edges file is an excel file rather than a csv file.
	if path.suffix in ['.xlsx', '.xls']:
		table = pandas.read_excel(path)
	elif path.suffix == '.csv':
		table = pandas.read_csv(path)
	elif path.suffix == '.tsv':
		table = pandas.read_table(path)
	else:
		message = f"The given filetype is not supported. Expected one of ['.csv', '.tsv', '.xls', '.xlsx'], received '{path.suffix}'"
		raise ValueError(message)
	return table

def main(edgesio:Union[Path, pandas.Series], annotations:Dict[str,str] = None, palette:Dict[str,str] = None, filename = None):
	"""
		Generates a lineageplot of the given data.
		Parameters
		----------
		edgesio: Union[Path, pandas.Series]
			The lineage for this workflow. If `edgesio` is a file it may contain additional information via these columns:
			- 'identity': The name of the genotype
			- `parent`: The name of the genotype's most curent ancestor
			- `annotation': An annotation to include alongside the genotype
			- `color`: The color to use when drawing the lineageplot
		annotations: Dict[str,str]
			An optional dictionary mapping annotations to the genotypes in `edges`. Overrides any annotations found in the `edges` file.
		palette: Dict[str,str]
			A colorscheme mapping genotype names to colors. May be given as a file, or extracted from the `edges` file.
		filename: Path
			The file to save the lineageplot as.
	"""

	if isinstance(edgesio, Path):
		table_edges = read_edges_filename(edgesio)
	else:
		# Assume `edgesio` is a pandas DataFrame
		table_edges = edgesio

	# Convert the table columns to lowercase so we don't have to work with multiple capitalizations.
	# Also convert to string in case other type managed to becom column labels.
	table_edges.columns = [str(i).lower for i in table_edges.columns]
	# Try to get the option from the input DataFrame.
	if annotations is None and 'annotation' in table_edges:
		annotations = table_edges['annotation'].to_dict()

	if palette is None:
		if 'color' in table_edges.column:
			# Use the given color.
			palette = table_edges['color'].to_dict()
		else:
			# Generate the palette
			pass


if __name__ == "__main__":
	pass