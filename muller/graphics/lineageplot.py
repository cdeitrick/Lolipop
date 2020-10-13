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


def flowchart(edges: pandas.DataFrame, palette: Dict[str, str], annotations: Dict[str, List[str]] = None, filename: Optional[Path] = None, add_score:bool = True)->pygraphviz.AGraph:
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
	identity_column = "Identity"
	parent_column = "Parent"

	# Check if the edges table is using the identity column an an index.
	if edges.index.name and edges.index.name == identity_column.lower():
		edges = edges.reset_index()
	if annotations is None:
		annotations = {}

	graph = pygraphviz.AGraph(splines = 'ortho', orientation = 90)
	graph.node_attr['shape'] = 'box'
	graph.node_attr['style'] = 'filled,rounded'
	graph.node_attr['fontname'] = 'lato'
	graph.edge_attr['dir'] = 'forward'

	for identity in edges[identity_column].unique():
		genotype_color = palette.get(identity)
		graph.add_node(identity, color = "#333333", fillcolor = genotype_color)
		node = graph.get_node(identity)
		node_label_properties = get_node_label_properties(identity, genotype_color, annotations.get(identity, []))
		node.attr.update(node_label_properties)

	for index, row in edges.iterrows():
		parent = row[parent_column]
		identity = row[identity_column]
		arguments = {
			'tooltip': f"{parent}",
			'labeldistance': 2.0
		}
		if add_score:
			score = row.get('score', 0)
			arguments['headlabel'] = f"{score:.1f}"
		graph.add_edge(parent, identity, **arguments)

	if filename:
		try:
			graph.draw(str(filename), prog = 'dot')
		except Exception as exception:
			logger.error(f"Could not generate lineage plot: {exception}")
	return graph
