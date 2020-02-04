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


def flowchart(edges: pandas.DataFrame, palette: Dict[str, str], annotations: Dict[str, List[str]] = None, filename: Optional[Path] = None)->pygraphviz.AGraph:
	"""
		Creates a lineage plot showing the ancestry of each genotype.
	Parameters
	----------
	edges: pandas.DataFrame
		Should have three columns: 'parent', 'identity', 'score'
	palette
	annotations
	filename: Optional[Path]
		Where to save the plot.

	Returns
	-------

	"""
	if annotations is None:
		annotations = {}

	graph = pygraphviz.AGraph(splines = 'ortho', orientation = 90)
	graph.node_attr['shape'] = 'box'
	graph.node_attr['style'] = 'filled,rounded'
	graph.node_attr['fontname'] = 'lato'
	graph.edge_attr['dir'] = 'forward'

	for identity in edges['Identity'].unique():
		genotype_color = palette.get(identity)
		graph.add_node(identity, color = "#333333", fillcolor = genotype_color)
		node = graph.get_node(identity)
		node_label_properties = get_node_label_properties(identity, genotype_color, annotations.get(identity, []))
		node.attr.update(node_label_properties)

	for index, row in edges.iterrows():
		parent = row['Parent']
		identity = row['Identity']
		score = row.get('score', 0)
		graph.add_edge(parent, identity, tooltip = f"Parent", headlabel = f"{score:.1f}", labeldistance = 2.0)

	if filename:
		try:
			graph.draw(str(filename), prog = 'dot')
		except Exception as exception:
			logger.error(f"Could not generate lineage plot: {exception}")
	return graph
