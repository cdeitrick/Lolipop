"""
Generates a palette using the information.
"""
from typing import Dict, List
import pandas

try:
	import colorset
	import treetools
except ModuleNotFoundError:
	import treetools
	from . import colorset


def apply_clade_colorscheme(clade: List[str], colorscheme: str) -> Dict[str, str]:
	"""
		Assigns colorschemes to each clade group
	Parameters
	----------
	clade: List[str]: A list of genotypes in this clade
	colorscheme: The name of a valid colorscheme

	Returns
	-------

	"""
	colorscheme = colorset.load_colorscheme(colorscheme, len(clade))
	palette = dict(zip(clade, colorscheme))
	return palette


def generate_lineage_palette(edges: pandas.DataFrame) -> Dict[str, str]:
	clades = treetools.parse_tree(edges)
	groups = clades.groupby(by = 'clade')

	palette = dict()
	for clade_label, colorscheme in zip(clades['clade'].unique(), colorset.distinctive_colorschemes):
		group = groups.get_group(clade_label)
		clade_palette = apply_clade_colorscheme(group.index, colorscheme)
		palette.update(clade_palette)
	return palette


if __name__ == "__main__":
	pass
