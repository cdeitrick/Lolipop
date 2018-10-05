from pathlib import Path
import numpy as np
import pandas as pd
from generate_muller_plot import generate_muller_dataframe, generate_muller_series
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import brewer

N = 20
cats = 10
df = pd.DataFrame(np.random.randint(10, 100, size=(N, cats))).add_prefix('y')

if __name__ == "__main__":

	folder = Path.home() / "Documents" / "github" / "muller_diagrams" / "Data files" / "B1_muller_try1"
	population = folder / "B1_muller_try1.ggmuller.populations.tsv"
	edges = folder / "B1_muller_try1.ggmuller.edges.tsv"
	script_filename = Path(__file__).with_name("r_script.r")
	output_filename = Path(__file__).with_name("df.csv")

	muller_df = generate_muller_dataframe(population, edges, script_filename, output_filename = output_filename)

	color_palette = [ '#333333',
		'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
		'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
		'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
		'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
		'#ffffff', '#000000'
	]
	colormap = dict()
	for color, genotype_label in zip(color_palette, sorted(muller_df['Identity'].unique())):
		colormap[genotype_label] = color

	muller_df['color'] = [colormap[i] for i in muller_df['Identity'].values]


	x, y, colors, labels = generate_muller_series(muller_df)

	areas = df
	x2 = np.hstack((df.index[::-1], df.index))

	p = figure(x_range = (0, 90), y_range = (0, 1))
	p.grid.minor_grid_line_color = '#eeeeee'


	print(y)
	x = [x]*len(y)
	print(x)
	p.patches([x] * len(y), y,
		color = colors, alpha = 0.8, line_color = None)

	output_file('brewer.html', title = 'brewer.py example')

	show(p)
