from pathlib import Path
from typing import Dict, List, Optional

import pandas
import pygraphviz
from loguru import logger
from muller import widgets


def get_node_label_properties(identity: str, genotype_color: str, annotation: List[str]) -> Dict[str, str]:
	luminance = widgets.calculate_luminance(genotype_color)
	if luminance > 0.5:
		font_color = "#333333"
	else:
		font_color = "#FFFFFF"

	label = [identity] + annotation
	label = "\n".join(label)

	label_properties = {
		'label':     label,
		'fontcolor': font_color
	}
	return label_properties


def flowchart(edges: pandas.DataFrame, palette: Dict[str, str], annotations: Dict[str, List[str]] = None, filenames: Optional[Dict[str,List[Path]]] = None):
	"""
		Creates a lineage plot showing the ancestry of each genotype.
	Parameters
	----------
	edges: pandas.DataFrame
		Should have three columns: 'parent', 'identity', 'score'
	palette
	annotations
	filenames: Optional[List[Path]]
		A list of filenames to save the flowchart as. Mainly useful if we want the output in multiple formats (ex. .png and .svg)

	Returns
	-------

	"""
	if annotations is None:
		annotations = {}
	if isinstance(filenames, Path):
		filenames = [filenames]

	graph = pygraphviz.AGraph(splines = 'ortho')
	graph.node_attr['shape'] = 'box'
	graph.node_attr['style'] = 'filled,rounded'
	graph.node_attr['fontname'] = 'lato'
	graph.edge_attr['dir'] = 'forward'

	for identity in edges['identity'].unique():
		genotype_color = palette.get(identity)
		graph.add_node(identity, color = "#333333", fillcolor = genotype_color)
		node = graph.get_node(identity)
		node_label_properties = get_node_label_properties(identity, genotype_color, annotations.get(identity, []))
		node.attr.update(node_label_properties)

	for index, row in edges.iterrows():
		parent = row['parent']
		identity = row['identity']
		score = row['score']
		graph.add_edge(parent, identity, tooltip = f"parent", headlabel = f"{score:.1f}", labeldistance = 2.0)
	if filenames:
		for filename in filenames:
			try:
				graph.draw(str(filename), prog = 'dot')
			except Exception as exception:
				logger.error(f"Could not generate lineage plot: {exception}")
	return graph
