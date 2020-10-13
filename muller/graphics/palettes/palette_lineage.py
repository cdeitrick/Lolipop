"""
Generates a palette using the information.
"""
import random
from typing import Dict, List

import pandas

try:
	from . import colorset
	from muller import treetools
except ModuleNotFoundError:
	from muller import treetools
	from muller.graphics.palettes import colorset


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


def generate_lineage_palette(edges: pandas.Series) -> Dict[str, str]:
	clades = treetools.parse_tree(edges)
	major_clades = get_major_clades(clades)
	clades['majorClade'] = [major_clades[i] for i in clades.index]
	groups = clades.groupby(by = 'majorClade')

	palette = dict()
	_schemes = colorset.distinctive_colorschemes + colorset.distinctive_colorschemes + colorset.distinctive_colorschemes # In case there are a very large number of clades.
	for clade_label, colorscheme in zip(clades['majorClade'].unique(), _schemes):
		if colorscheme is None:
			colorscheme = random.choice(colorset.distinctive_colorschemes)
		group = groups.get_group(clade_label)
		clade_palette = apply_clade_colorscheme(group.index, colorscheme)
		palette.update(clade_palette)
	return palette


def get_major_clades(tree_table: pandas.DataFrame):
	clade_counts = tree_table['Parent'].value_counts()
	cutoff = int(len(tree_table) / 20)
	major_clades = list(clade_counts[clade_counts > cutoff].index)
	major_clade_map = dict()
	for child in tree_table.index:
		if child in major_clades:
			first_major_clade = child
		else:
			parents = treetools.get_parent_nodes(tree_table, child)
			try:
				first_major_clade = [i for i in parents if i in major_clades][0]
			except IndexError:
				first_major_clade = child
		major_clade_map[child] = first_major_clade
	return major_clade_map
