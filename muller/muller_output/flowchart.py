import pygraphviz
import dataio
import pandas
import widgets
from pathlib import Path
from typing import Dict, List, Any
def flowchart(edges:pandas.DataFrame, palette:Dict[str,str], clusters:Any = None, annotations:Dict[str,List[str]] = None, filename:Path = None):
	graph = pygraphviz.AGraph(splines = 'ortho')
	graph.node_attr['shape'] = 'box'
	graph.node_attr['style'] = 'filled,rounded'
	graph.node_attr['fontname'] = 'lato'
	graph.edge_attr['dir'] = 'forward'

	for identity in edges['Identity'].unique():
		genotype_color = palette.get(identity)
		graph.add_node(identity, color = "#333333", fillcolor = genotype_color)

		if annotations:
			luminance = widgets.calculate_luminance(genotype_color)
			if luminance > 0.5:
				font_color = "#333333"
			else:
				font_color = "#FFFFFF"
			node = graph.get_node(identity)
			label = [identity] + annotations.get(identity, [])
			label = "\n".join(label)

			node.attr['label'] = label
			node.attr['fontcolor'] = font_color



	for index, row in edges.iterrows():
		parent = row['Parent']
		identity = row['Identity']

		graph.add_edge(parent,identity, tooltip = f"parent")
	graph.draw(str(filename), prog = 'dot')

if __name__ == "__main__":
	pass

