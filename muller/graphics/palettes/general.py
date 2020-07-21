from typing import *

import pandas

from . import colorset
from .palette_annotation import generate_annotation_palette
from .palette_distinctive import generate_distinctive_palette
from .palette_lineage import generate_lineage_palette


def generate_palette(edges: pandas.Series, custom_palette: Dict[str, str] = None, annotations: Dict[str, List[str]] = None,
		kind = 'distinctive') -> Dict[str, str]:
	"""
	Parameters
	----------
	edges: A list of unique genotypes or the edges table.
	custom_palette: Dict[str,str]
	annotations
	kind
	"""
	if custom_palette is None: custom_palette = {}
	if annotations is None: custom_palette = {}

	unique_genotypes = list(edges.index)

	if kind == 'lineage':
		palette = generate_lineage_palette(edges)
	elif kind == 'annotation':
		palette = generate_annotation_palette(annotations, custom_palette)
	else:
		palette = generate_distinctive_palette(unique_genotypes)
	for label in unique_genotypes:
		new_color = custom_palette.get(label, palette.get(label))
		if new_color is None:
			palette[label] = colorset.random_color()

	palette['genotype-0'] = "#FFFFFF"
	palette['removed'] = "#333333"
	return palette
