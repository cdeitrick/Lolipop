from muller import dataio
from muller.graphics import lineageplot


def test_get_font_properties():
	expected = {'label': "genotype-1\ngene1\ngene2", 'fontcolor': '#FFFFFF'}
	assert lineageplot.get_node_label_properties('genotype-1', '#330033', ['gene1', 'gene2']) == expected

	expected = {'label': 'genotype-11', 'fontcolor': '#333333'}
	assert lineageplot.get_node_label_properties('genotype-11', '#FFFFFF', []) == expected


def test_flowchart():
	edges_string = """
	Parent	Identity	score
	genotype-0	genotype-1	1
	genotype-1	genotype-13	2
	"""
	edges_table = dataio.import_table(edges_string)
	palette = {"genotype-13": '#222222', 'genotype-1': '#CCCCCC', 'genotype-0': '#000000'}

	resultgraph = lineageplot.flowchart(edges_table, palette, annotations = {'genotype-1': ['gene1']})

	dark_node = resultgraph.get_node('genotype-13')
	assert dark_node.attr['fontcolor'] == '#FFFFFF'
	assert dark_node.attr['label'] == 'genotype-13'

	light_node = resultgraph.get_node('genotype-1')
	assert light_node.attr['fontcolor'] == '#333333'
	assert light_node.attr['label'] == 'genotype-1\ngene1'
