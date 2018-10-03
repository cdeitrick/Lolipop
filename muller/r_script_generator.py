from matplotlib import pyplot as plt
import pandas
from pathlib import Path
import subprocess
from typing import Dict

plt.style.use('seaborn')
import altair
from plotnine import *

pandas.set_option('display.max_rows', 500)
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 250)


def execute_r_script(path: Path, script: str) -> Path:
	path.write_text(script)
	print("generating muller plot...")
	process = subprocess.run(
		['Rscript', '--vanilla', '--silent', path],
		stdout = subprocess.PIPE,
		stderr = subprocess.PIPE
	)
	print(process.stdout)
	print(process.stderr)
	_extra_file = Path.cwd() / "Rplots.pdf"
	if _extra_file.exists():
		_extra_file.unlink()
	return path


def generate_muller_dataframe(population: Path, edges: Path, script_filename: Path, output_filename: Path) -> pandas.DataFrame:
	script = """
		library(ggplot2)
		library(ggmuller)

		population <- read.table("{population}", header=TRUE)
		edges <- read.table("{edges}", header=TRUE)

		Muller_df <- get_Muller_df(edges, population)

		write.csv(Muller_df, "{path}")
		""".format(
		path = output_filename.absolute(),
		population = population.absolute(),
		edges = edges.absolute(),
	)

	execute_r_script(script_filename, script)
	# if script_filename.exists():
	#	script_filename.unlink()

	muller_df = pandas.read_csv(output_filename)
	return muller_df


def generate_r_script(annotations: Dict[str, str], trajectory: Path, population: Path, edges: Path, output_file: Path,
		color_palette: Dict[str, str]) -> str:
	an = ["-".join(v) for k, v in sorted(annotations.items())]
	an = ['"{}"'.format(i) for i in an]
	script = """
	library(ggplot2)
	library(ggmuller)

	population <- read.table("{population}", header=TRUE)
	edges <- read.table("{edges}", header=TRUE)
	trajectories <- read.table("{trajectory}", header = TRUE, sep = "\t")

	Muller_df <- get_Muller_df(edges, population)
	#Muller_plot(Muller_df, add_legend = TRUE)

	palette <- c({palette})

	text <- c({labels})

	ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
    geom_area() +
    geom_text() + 
    theme(legend.position = "right") +
    guides(linetype = FALSE, color = FALSE) + 
    scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_fill_manual(name = "Identity", values = palette) +
    scale_color_manual(values = palette) +

    annotate(c, x = 2:5, y = 25, label = "Some text")

	ggsave("{output}", height = 10, width = 10)
	""".format(
		labels = ",".join(an),
		trajectory = trajectory,
		population = population.absolute(),
		edges = edges.absolute(),
		output = output_file.absolute(),
		palette = ",".join(['"#333333"'] + ['"{}"'.format(v) for k, v in sorted(color_palette.items())])
	)
	script = '\n'.join(i.strip() for i in script.split('\n'))
	print(script)
	return script


def generate_muller_plot(muller_df: pandas.DataFrame):
	color_palette = [
		'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
		'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
		'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
		'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
		'#ffffff', '#000000'
	]
	length = muller_df['Identity'].nunique()
	palette = color_palette[: length]
	fig, ax = plt.subplots(figsize = (10, 10))
	print(muller_df)

	gplot = (
			ggplot(muller_df, aes(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", color = "Identity")) +
			geom_area() +
			theme(legend_position = "right") +
			# guides(scale_linetype = False, color = False)+
			scale_y_continuous(name = "Percentage", expand = (0, 0)) +
			scale_x_continuous(name = "Generation", expand = (0, 0)) +
			scale_fill_manual(name = "Identity", values = palette) +
			scale_color_manual(values = palette)
	)
	print(gplot)


def generate_new_muller_plot(muller_df: pandas.DataFrame):
	color_palette = [
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

	length = muller_df['Identity'].nunique()
	palette = color_palette[: length]
	fig, ax = plt.subplots(figsize = (10, 10))
	genotype_order = list()
	for index, row in muller_df.iterrows():
		genotype_label = row['Group_id']
		if genotype_label not in genotype_order:
			genotype_order.append(genotype_label)
	gens = list()
	x = list()
	y = list()
	colors = list()
	for genotype_label in genotype_order:
		ddf = muller_df[muller_df['Group_id'] == genotype_label]
		color = ddf['color'].tolist()[0]
		colors.append(color)
		x = ddf['Generation'].tolist()
		y.append(ddf['Frequency'].tolist())

	plt.stackplot(x, y, colors = colors)

	# Basic stacked area chart.

	plt.legend(loc = 'upper left')
	plt.show()


def generate_altair(muller_df: pandas.DataFrame):
	color_palette = [
		'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
		'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
		'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
		'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
		'#ffffff', '#000000'
	]
	colormap = dict()
	for color, genotype_label in zip(color_palette, muller_df['Identity'].unique()):
		colormap[genotype_label] = color

	muller_df['color'] = [colormap[i] for i in muller_df['Identity'].values]

	length = muller_df['Identity'].nunique()
	palette = color_palette[: length]
	fig, ax = plt.subplots(figsize = (10, 10))

	chart = altair.Chart(muller_df).mark_area().encode(
		x = 'Generation',
		y = 'Frequency',
		color = 'color'
	)
	chart.save('chart.html')

	plt.show()


if __name__ == "__main__":
	folder = Path.home() / "Documents" / "github" / "muller_diagrams" / "Data files" / "B1_muller_try1"
	population = folder / "B1_muller_try1.ggmuller.populations.tsv"
	edges = folder / "B1_muller_try1.ggmuller.edges.tsv"
	script_filename = Path(__file__).with_name("r_script.r")
	output_filename = Path(__file__).with_name("df.csv")

	df = generate_muller_dataframe(population, edges, script_filename, output_filename = output_filename)
	# df = df.set_index('Unnamed: 0')
	# df = df.sort_values(by = ['Generation', 'Unique_id'], ascending = True)

	generate_new_muller_plot(df)
