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
		'M1/M2': ["yciR Y355D", "rpfR E481D"],
		'M3': ["OGHD R204S"],
		'M4/M5': ["LysR-like L40V", "d95 d95"],
		'M5': ["d95 d95"],
		'M6': ["manC STOP"],
		'M7': ["bfr 5`"],
		'M8': ["succinate dehydrogenase synonymous"],
		'M9': ["d49 d49"],
		'M10':["DUF88 A209P"],
		'M11': ["MltA 5`"],
		'M12': ["wspD L35P"],
		'M13': ["wspA A407V"],
		'M14': ["wspE S726L"],
		'M15': ["rpfR A106P"],
		'M16': ["rpoC A1318V"],
		'M17': ["wspA I196N"],
		'M18': ["bfr 5`"],
		'M19': ["McsS Synonymous"],
		'M20': ["manC T154A"]
	}
	other = "#111111"
	genotypes = [
		"#fca082",
		'#37a055',
		'#9e9ac8',
		'#ff5c00',
		'#3787c0',
		'#ff7f00'
	]
	singleton = "#FFFFFF"
	palette = {
		'M1/M2': genotypes[0],
		'M3': '#aedea7',
		'M4/M5': genotypes[1],
		'M6': genotypes[2],
		'M7': genotypes[2],
		'M8': genotypes[3],
		'M9': '#ffff07',
		'M10': genotypes[3],
		'M11': genotypes[3],
		'M12': genotypes[3],
		'M13': other,
		'M14': other,
		'M15': genotypes[4],
		'M16': genotypes[4],
		'M17': genotypes[4],
		'M18': other
	}
	output = Path("/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/muller/manuscript/figures/results/flowchart.png")

	flowchart(edges, annotations = annotations, palette = palette, filename = output)