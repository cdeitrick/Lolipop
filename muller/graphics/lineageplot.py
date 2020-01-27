from pathlib import Path
from typing import Dict, List, Optional

import pandas
import pygraphviz
from loguru import logger
from muller import widgets
from io import StringIO

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

if __name__ == "__main__":
	# Identity : Parent
	edges = {
		'M1/M2':  'genotype-0',
		'M3':  "M1/M2",
		'M4/M5': "M3",
		'M6':  "M4/M5",
		'M7':  "M6",
		'M8':  "M7",
		'M9':  "M8",
		'M10': "M7",
		'M11': "M10",
		'M12': "M11",
		'M13': "M7",
		'M14': "M7",
		'M15': "genotype-0",
		'M16': "genotype-0",
		'M17': "M16",
		'M18': "M17"
		#'M19': ["McsS Synonymous"],
		#'M20': ["manC T154A"]
	}
	identity, parent = zip(*sorted(edges.items()))
	edges = pandas.DataFrame(
		{
			'Identity': identity,
			'Parent': parent
		}
	)
	annotations = {
		'M1/M2': ["yciR Y355D", "monoxygenase E481D"],
		'M3': ["OGHD R204S"],
		'M4/M5': ["LysR-like", "95 gene deletion"],
		'M6': ["manC -1bp"],
		'M7': ["bfr 5` -10bp"],
		'M8': ["49 gene deletion"],
		'M9': ["succinate dehydrogenase G147G"],
		'M10':["DUF88 A209P"],
		'M11': ["MltA 5` -10bp"],
		'M12': ["wspD L35P"],
		'M13': ["wspA A407V"],
		'M14': ["wspE S726L"],
		'M15': ["yciR A106P"],
		'M16': ["wspA I196N"],
		'M17': ["rpoC A1318V"],
		'M18': ["bfr 5` -35bp"],
		'M19': ["McsS R135R"],
		'M20': ["manC A412V"]
	}

	genotypemap = [
		["M1/M2"],
		["M13", "M20"],
		["M15", "M16", "MM17"],
		["M14", "M18"],
		["M3"],
		["M4/M5"],
		["M6", "M7"],
		["M19"],
		["M10", "M11", "M12", "M9"],
		["M8"]
	]
	singleton = "#FFFFFF"
	palette = {
		'M1/M2': "#fca082",
		'M3': "#aedea7",
		'M4/M5': "#371055",
		'M6': "#9e9ac8",
		'M7': "#9e9ac8",
		'M8': "#ffff07",
		'M9': "#ff5c00",
		'M10': "#ff5c00",
		'M11': "#ff5c00",
		'M12': "#ff5c00",
		'M13': "#adb0e6",
		'M14': "#e32f27",
		'M15': "#3787c0",
		'M16': "#3787c0",
		'M17': "#3787c0",
		'M18': "#e32f27",
		'M20': "#adb0e6"
	}
	output = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/muller/manuscript/figures/results/flowchart.png")

	flowchart(edges, annotations = annotations, palette = palette, filename = output)