import subprocess
from pathlib import Path
from typing import Dict, List


def generate_r_script(trajectory: Path, population: Path, edges: Path, plot_filename: Path,
		color_palette: Dict[str, str], genotype_labels: List[str]) -> str:
	# `color_palette` should be an OrderedDict.
	genotype_labels = sorted(genotype_labels)

	script_colors = ",".join(['"{}"'.format(color_palette.get(k, "#333333")) for k in genotype_labels])

	script = """
	library(ggplot2)
	library(ggmuller)

	population <- read.table("{population}", header=TRUE)
	edges <- read.table("{edges}", header=TRUE)

	Muller_df <- get_Muller_df(edges, population)
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
		palette = script_colors
	)
	script = '\n'.join(i.strip() for i in script.split('\n'))

	return script


def execute_r_script(path: Path, script: str) -> Path:
	""" Executes `script` by first saving it to `path` then calling `Rscript` to run it."""
	path.write_text(script)
	subprocess.run(
		['Rscript', '--vanilla', '--silent', path],
		stdout = subprocess.PIPE,
		stderr = subprocess.PIPE
	)
	# Running the script will generate an Rplots.pdf file which we don't care about.
	# `ggsave` is used to save the plot instead.
	_extra_file = Path.cwd() / "Rplots.pdf"
	if _extra_file.exists():
		_extra_file.unlink()
	return path
