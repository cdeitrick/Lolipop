import subprocess
from pathlib import Path
from typing import Dict, Optional, List

import pandas


def generate_mermaid_diagram(backgrounds: pandas.DataFrame, color_palette: Dict[str, str]) -> str:
	"""
	graph LR
    id1(Start)-->id2(Stop)
    style id1 fill:#f9f,stroke:#333,stroke-width:4px
    style id2 fill:#ccf,stroke:#f66,stroke-width:2px,stroke-dasharray: 5, 5
	Parameters
	----------
	backgrounds
	color_palette

	Returns
	-------

	"""
	genotype_labels = list(set(backgrounds['Identity'].values))

	node_map = {k: "style id{} fill:{}".format(k.split('-')[-1], color_palette[k]) for k in genotype_labels}

	diagram_contents = ["graph TD;"]
	for _, row in backgrounds.iterrows():
		parent = row['Parent']
		identity = row['Identity']
		parent_id = parent.split('-')[-1]
		identity_id = identity.split('-')[-1]
		line = "id{left_id}({left})-->id{right_id}({right});".format(
			left_id = identity_id,
			right_id = parent_id,
			left = identity,
			right = parent
		)
		diagram_contents.append(line)

	diagram_contents += list(node_map.values())

	script_text = "\n".join(diagram_contents)
	return script_text


def generate_r_script(trajectory: Path, population: Path, edges: Path, table_filename: Path, plot_filename: Path, script_filename: Path,
		color_palette: Dict[str, str], genotype_labels:List[str]) -> Optional[pandas.DataFrame]:

	# `color_palette` should be an OrderedDict.
	genotype_labels = sorted(genotype_labels)
	script_colors = ",".join(['"{}"'.format(color_palette[k]) for k in genotype_labels])

	script = """
	library(ggplot2)
	library(ggmuller)

	population <- read.table("{population}", header=TRUE)
	edges <- read.table("{edges}", header=TRUE)

	Muller_df <- get_Muller_df(edges, population)
	write.csv(Muller_df, "{path}", sep = "\\t", col.names = NA)
	palette <- c({palette})

	ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
    geom_area() +
    theme(legend.position = "right") +
    guides(linetype = FALSE, color = FALSE) + 
    scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_fill_manual(name = "Identity", values = palette) +
    scale_color_manual(values = palette)

	ggsave("{output}", height = 10, width = 10)
	""".format(
		trajectory = trajectory,
		population = population.absolute(),
		edges = edges.absolute(),
		output = plot_filename.absolute(),
		palette = script_colors,
		path = table_filename
	)
	script = '\n'.join(i.strip() for i in script.split('\n'))

	execute_r_script(script_filename, script)
	if table_filename.exists():
		muller_df = pandas.read_csv(table_filename)
	else:
		muller_df = None
	return muller_df


def execute_r_script(path: Path, script: str) -> Path:
	path.write_text(script)
	print("generating muller plot...")
	subprocess.run(
		['Rscript', '--vanilla', '--silent', path],
		stdout = subprocess.PIPE,
		stderr = subprocess.PIPE
	)
	_extra_file = Path.cwd() / "Rplots.pdf"
	if _extra_file.exists():
		_extra_file.unlink()
	return path


def excecute_mermaid_script(path: Path, script:str, mermaid_render: Path):
	path.write_text(script)
	try:
		subprocess.call(
			["mmdc", "--height", "400", "-i", path, "-o", mermaid_render],
			stdout = subprocess.PIPE,
			stderr = subprocess.PIPE)
	except FileNotFoundError:
		pass
