#!python
"""
	Main script to run the muller workflow.
	.
	|---- .clade.annotated.(png|svg)
	|---- .genotypes.tsv
	|---- .lineage.clade.(png|svg)
	|---- figures/
	|----|---- clade/
	|----|----|---- .clade.annotated.(png|svg)
	|----|----|---- .clade.unannotated.(png|svg)
	|----|----|---- .clade.genotypes.(png|svg)
	|----|----|---- .clade.trajectories.(png|svg)
	|----|---- unique/
	|----|----|---- .unique.annotated.(png|svg)
	|----|----|---- .unique.unannotated.(png|svg)
	|----|----|---- .unique.genotypes.(png|svg)
	|----|----|---- .unique.trajectories.(png|svg)
	|----|---- .dendrogram.png
	|----|---- .pairwisedistance.svg
	|---- scripts/
	|----|---- .r
	|---- suplementary-files/
	|----|---- .genotypeinformation.json
	|----|---- .lineagescores.tsv
	|----|---- .options.json
	|---- tables/
	|----|---- .trajectories.tsv
	|----|---- .genotypes.tsv
	|----|---- .distance.tsv
	|----|---- .linkagematrix
	|----|---- .ggmuller.edges.tsv
	|----|---- .ggmuller.populations.tsv
	|----|---- .mullerdataframe.tsv
"""
import argparse
from pathlib import Path
from typing import *

import pandas

from muller import graphics, widgets
from muller.graphics import graphicsio

pandas.set_option('mode.chained_assignment',
	None)  # This disables the warning about setting a value on a copy of a dataframe.
from loguru import logger

from muller import clustering, dataio, inheritance, commandline_parser
from muller.dataio import projectdata, annotations, projectpaths

logger.remove()  # Need to remove the default sink so that the logger doesn't print messages twice.
import sys

if commandline_parser.DEBUG:
	logger.add(sys.stderr, level = "DEBUG")
else:
	logger.add(sys.stderr, level = 'INFO', format = "{time:YYYY-MM-DD HH:mm:ss} {level} {message}")


def run_genotype_inference_workflow(trajectoryio: Union[str, Path, pandas.DataFrame], metric: str, dlimit: float,
		flimit: float,
		similarity_cutoff: float, known_genotypes: Optional[Path] = None, threads: Optional[int] = None,
		is_genotype: bool = False) -> projectdata.DataGenotypeInference:
	"""
	Parameters
	----------
	trajectoryio: Union[str,Path,pandas.DataFrame]
		Either the path to a table with mutational trajectories or the table itself.
		The table should already be cleaned/filtered.
		The table should be organized with timepoints as columns and rows as genotype frequencies.
		Ex.
		genotype	0	7	13 ...
		genotype-1	0	.1	.3 ...
		genotype-2	0.2	.1 ...
	metric: Literal['']
	dlimit, flimit: float
	similarity_cutoff: float
	known_genotypes
	threads
	is_genotype: bool
	"""
	if isinstance(trajectoryio, (str, Path)):
		logger.info(f"Reading '{trajectoryio}' as the trajectory table.")
		trajectories, trajectory_info = dataio.parse_trajectory_table(trajectoryio)
	else:
		trajectories = trajectoryio
		trajectory_info = None

	# Read in the starting genotypes.
	if known_genotypes:
		known_genotypes = dataio.parse_known_genotypes(known_genotypes)

	genotype_generator = clustering.ClusterMutations(
		metric = metric,
		dlimit = dlimit,
		flimit = flimit,
		starting_genotypes = known_genotypes,
		threads = threads
	)
	if is_genotype:
		logger.info(f"Skipping genotype infeerence...")
		table_genotypes = trajectories.copy(deep = True)
		table_genotypes = genotype_generator.organizer.run(table_genotypes)
		# Make sure the "Genotype" column is properly named.
		table_genotypes.index.name = "Genotype"
		genotype_members = {f"genotype-{k}": [k] for k in table_genotypes.index}
		# Make sure the genotype names are prefixed with 'genotype'
		table_genotypes.index = [f'genotype-{i}' for i in table_genotypes.index]
		genotype_data = projectdata.DataGenotypeInference(
			table_trajectories = trajectories,
			table_genotypes = table_genotypes,
			genotype_members = genotype_members,
			matrix_distance = None,
			clusterdata = None,
			table_trajectories_info = None
		)
	else:
		genotype_data = genotype_generator.run(trajectories, distance_cutoff = similarity_cutoff)
	genotype_data.table_trajectories_info = trajectory_info
	return genotype_data


def run_genotype_lineage_workflow(genotypeio: Union[str, Path, pandas.DataFrame], dlimit: float, flimit: float,
		pvalue: float, known_ancestry: Optional[Path], conservative:bool) -> projectdata.DataGenotypeLineage:
	"""

	Parameters
	----------
	genotypeio
	dlimit
	flimit
	pvalue
	known_ancestry

	Returns
	-------

	"""
	lineage_generator = inheritance.LineageWorkflow(
		dlimit = dlimit,
		flimit = flimit,
		pvalue = pvalue,
		conservative = conservative
	)

	# Read in the input data if it is not already a pandas.DataFrame object
	if isinstance(genotypeio, (str, Path)):
		logger.info(f"Reading '{genotypeio}' as the genotype table.")
		genotypes = dataio.parse_genotype_table(genotypeio)
	else:
		genotypes = genotypeio

	if known_ancestry:
		known_ancestry = dataio.read_map(known_ancestry)
		# Remove the table headers.
		if 'Parent' in known_ancestry:
			known_ancestry.pop('Parent')
		if 'Identity' in known_ancestry:
			known_ancestry.pop('Identity')

	lineage_data = lineage_generator.run(genotypes, known_ancestry)
	return lineage_data


def run_lineageplot_workflow(edgesio: Union[str, Path, pandas.DataFrame], filename: Path,
		sheet_name: Optional[str] = None):
	""" Generates a lineage plot given an `edges` table. The columns should be named `parent` and `identity`.
		An optional `annotation` column will be used to annotate the plot.
	"""
	from muller.graphics import flowchart
	from muller.graphics.palettes import generate_palette
	if isinstance(edgesio, (str, Path)):
		edges = dataio.import_table(edgesio, sheet_name = sheet_name)
	else:
		edges = edgesio
	# Convert the column labels to lowercase
	edges.columns = [i.capitalize() for i in edges.columns]
	edges = edges.set_index('Identity')
	if 'Annotation' in edges.columns:
		annotations = edges.pop("Annotation").to_dict()
		# Make sure 'annotations' is a list of str which is expected by the flowchart
		annotations = {k: [v] for k, v in annotations.items()}
	else:
		annotations = dict()

	palette = generate_palette(edges['Parent'], kind = 'lineage')

	lineage = flowchart(edges, palette, annotations, filename)


def get_base_filename(filename: Path, name: Optional[str], sheetname: Optional[str]) -> str:
	if name:
		base_filename = name
	else:
		base_filename = filename.stem
	if sheetname and sheetname != 'Sheet1':
		base_filename += '.' + sheetname

	return base_filename

class WorkFlowLineage:
	def __init__(self, program_options):
		self.program_options = program_options

	def run(self, filename:Path):
		pass
def run_workflow(program_options: argparse.Namespace):

	# TODO: Test whether the lineage makes sense by computing the sum of genotypes/lineages at each timepoint,
	# Where 100% should be the maximum value.

	"""
		TODO: Running this command will fail:
		python /home/cld100/Documents/github/muller_diagrams/lolipop lineage
		--input traverse-etal-B1-mutationfrequencies.xlsx
		--sheetname trajectories.filtered
		--output traverse-genotypes
		--known-ancestry known_ancestry2.txt
		--genotypes
	"""

	logger.info("Running with the following parameters")
	for key, value in vars(program_options).items():
		logger.info(f"\t{key:<30}{value}")
	output_folder = program_options.output_folder
	paths = projectpaths.OutputFilenames(output_folder, program_options.filename.stem)
	logger.info("Parsing options...")
	data_basic = projectdata.DataWorkflowBasic(
		version = commandline_parser.__VERSION__,
		filename = program_options.filename,
		program_options = program_options
	)
	if not data_basic.program_options.output_folder.exists():
		data_basic.program_options.output_folder.mkdir()
	logger.info("Importing trajectories...")
	trajectory_table, trajectory_info = dataio.parse_trajectory_table(program_options.filename,
		program_options.sheetname)

	# Need to read in the input dataset.
	logger.info("Generating genotypes...")
	result_genotype_inference = run_genotype_inference_workflow(
		trajectory_table,
		program_options.metric,
		dlimit = program_options.dlimit,
		flimit = program_options.flimit,
		similarity_cutoff = program_options.similarity_cutoff,
		known_genotypes = program_options.known_genotypes,
		threads = program_options.threads,
		is_genotype = program_options.is_genotype
	)

	if result_genotype_inference.table_trajectories_info is None:
		result_genotype_inference.table_trajectories_info = trajectory_info

	genotype_annotations = annotations.read_genotype_annotations(trajectory_info,
		result_genotype_inference.genotype_members)
	result_genotype_lineage = run_genotype_lineage_workflow(
		result_genotype_inference.table_genotypes,  # should already be sorted.
		dlimit = program_options.dlimit,
		flimit = program_options.flimit,
		pvalue = program_options.pvalue,
		known_ancestry = program_options.known_ancestry,
		conservative = program_options.conservative
	)

	paths.save_projectdata_basic(data_basic)
	paths.save_workflow_clustering(result_genotype_inference)
	paths.save_workflow_hierarchy(result_genotype_inference.clusterdata)
	paths.save_workflow_lineage(result_genotype_lineage)

	# save_tables(data_basic, result_genotype_inference, result_genotype_lineage, genotype_annotations)
	# Save using the older graphics workflow for now.
	render_graphics(
		paths = paths,
		data_basic = data_basic,
		data_inference = result_genotype_inference,
		data_lineage = result_genotype_lineage,
		genotype_annotations = genotype_annotations
	)

	data_basic.save(output_folder)


def save_tables(data_basic: projectdata.DataWorkflowBasic, data_inference: projectdata.DataGenotypeInference,
		data_lineage: projectdata.DataGenotypeLineage):
	"""
		.figures
		|---- tables
		|----|---- genotypes.tsv
		|----|---- trajectories.tsv
		|----|---- trajectores.filtered.tsv
		|----|---- edges.tsv
		|----|---- populations.tsv
		|----|---- distances.tsv
		|----|---- scorerecords.tsv
		|----|---- linkagetable.tsv
		|----|---- parameters.json
		|----|---- genotypeinformation.json
		|---- scripts
		|----|---- lineage.r
	"""
	delimiter = "\t"
	suffix = "tsv"
	prefix = get_base_filename(data_basic.program_options.filename, data_basic.program_options.name,
		data_basic.program_options.sheetname)
	folder_tables = widgets.checkdir(data_basic.program_options.output_folder / "tables")
	folder_data = widgets.checkdir(data_basic.program_options.output_folder / "data")
	data_basic.save(folder_data)
	data_inference.save(folder_tables, prefix)
	data_lineage.save(folder_tables, prefix)
	"""
	filename_table_genotypes = folder_table / f"{prefix}.genotypes.{suffix}"
	filename_table_trajectories = folder_table / f"{prefix}.trajectories.{suffix}"
	filename_table_trajectories_filtered = folder_table / f"{prefix}.trajectories.filtered.{suffix}"
	filename_table_edges = folder_table / f"{prefix}.edges.{suffix}"
	filename_table_population = folder_table / f"{prefix}.population.{suffix}"
	filename_table_distances = folder_table / f"{prefix}.distances.{suffix}"
	filename_table_scorerecords = folder_table / f"{prefix}.scorerecords.{suffix}"
	filename_table_linkage = folder_table / f"{prefix}.linkage.{suffix}"

	filename_data_parameters = folder_data / f"{prefix}.parameters.json"
	filename_data_genotypes = folder_data / f"{prefix}.genotypes.json"

	data_inference.table_genotypes.to_csv(filename_table_genotypes, sep = delimiter)
	data_inference.table_trajectories.to_csv(filename_table_trajectories, sep = delimiter)
	data_inference.matrix_distance.squareform().to_csv(filename_table_distances, sep = delimiter)
	data_inference.clusterdata.table_linkage.to_csv(filename_table_linkage, sep = delimiter)
	
	data_lineage.table_edges.to_csv(filename_table_edges, sep = delimiter)
	data_lineage.table_populations.to_csv(filename_table_population, sep = delimiter)
	data_lineage.table_scores.to_csv(filename_table_scorerecords, sep = delimiter)
	"""


def render_graphics(paths: projectpaths.OutputFilenames, data_basic: projectdata.DataWorkflowBasic,
		data_inference: projectdata.DataGenotypeInference,
		data_lineage: projectdata.DataGenotypeLineage, genotype_annotations: Dict[str, List[str]]):
	""" Graphics are parametrized by filename suffix (png vs svg) and palette.
		.figures
		|---- dendrogram.png
		|---- heatmap.png
		|---- distancedistribution.png
		|---- unique palette
		|    |---- muller.annotated.(svg|png)
		|    |---- muller.unannotated.(svg|png)
		|    |---- timeseriespanel.(svg|png)
		|    |---- lineageplot.(svg|png)
		|---- lineage palette
		|    |---- muller.annotated.(svg|png)
		|    |---- muller.unannotated.(svg|png)
		|    |---- timeseriespanel.(svg|png)
		|    |---- lineageplot.(svg|png)
	"""
	custom_palette = dataio.read_map(data_basic.program_options.genotype_palette_filename)
	prefix = get_base_filename(
		data_basic.program_options.filename,
		data_basic.program_options.name,
		data_basic.program_options.sheetname
	)
	workflow_graphics = graphicsio.OutputGeneratorGraphics(
		project_folder = paths.folder_output,
		palettes = [],
		sample_basename = prefix,
		render = True
	)

	# Plot the figures that aren't parametrized
	if data_inference.clusterdata is not None:
		workflow_graphics.generate_dendrogram(data_inference.clusterdata.table_linkage, data_inference.matrix_distance,
			paths.filename_figure_linkage_plot)
	if data_inference.matrix_distance is not None:
		workflow_graphics.generate_heatmap(data_inference.matrix_distance.squareform(),
			paths.filename_figure_distance_heatmap)
	if data_inference.clusterdata is not None and data_inference.matrix_distance is not None:
		graphics.generate_distance_plot(
			data_inference.matrix_distance.values,
			data_inference.clusterdata.distance_cutoff,
			paths.filename_figure_distribution
		)
	# Set up the generators
	generator_panel_timeseries = graphics.TimeseriesPanel(render = data_basic.program_options.render)
	generator_plot_timeseries = graphics.TimeseriesPlot(render = data_basic.program_options.render)
	generator_plot_muller = graphics.MullerPlot(outlines = data_basic.program_options.draw_outline,
		render = data_basic.program_options.render)
	generator_plot_muller_panel = graphics.MullerPanel()
	# These figures need to be parametrized by palette type and filetype.
	for palette_name in ["unique", "lineage"]:
		folder_palette = paths.folder_figures_lineage if palette_name == 'lineage' else paths.folder_figures_unique

		current_palette_data = graphics.palettes.generate_palette(
			edges = data_lineage.table_edges,
			custom_palette = custom_palette,
			annotations = genotype_annotations,
			kind = palette_name,
		)

		current_palette = graphics.Palette(
			name = palette_name,
			palette = current_palette_data,
			members = data_inference.genotype_members
		)

		filename_palette = paths.filename_palette.with_name(f"palette.{palette_name}.json")
		current_palette.save(filename_palette)

		for suffix in ["png", "svg"]:
			filename_mullerplot_annotated = paths.get_template(folder_palette,
				paths.template_figure_muller_diagram_annotated, suffix)
			filename_mullerplot_unannotated = paths.get_template(folder_palette,
				paths.template_figure_muller_diagram_unannotated, suffix)
			filename_timeseries_panel = paths.get_template(folder_palette, paths.template_figure_panel_timeseries,
				suffix)
			filename_timeseries_genotypes = paths.get_template(folder_palette,
				paths.template_figure_timeseries_genotype, suffix)
			# filename_mullerplot_annotated = folder_palette / f"{prefix}.muller.annotated.{suffix}"
			# filename_mullerplot_unannotated = folder_palette / f"{prefix}.muller.unannotated.{suffix}"
			# filename_timeseries_panel = folder_palette / f"{prefix}.timeseriespanel.{suffix}"
			# filename_timeseries_genotypes = folder_palette / f"{prefix}.timeseries.genotypes.{suffix}"
			# filename_lineageplot = folder_palette / f"{prefix}.lineageplot.{suffix}"
			# filename_muller_panel = folder_palette / f"{prefix}.mullerpanel.{suffix}"
			filename_graphviz_script = folder_palette / f"{prefix}.lineagescript.{suffix}.dot"
			filename_lineage_panel = folder_palette / f"{prefix}.lineagepanel.{suffix}"

			filename_lineageplot = paths.get_template(folder_palette, paths.template_figure_lineageplot, suffix)
			filename_muller_panel = paths.get_template(folder_palette, paths.template_figure_panel_muller, suffix)
			# filename_script_graphviz = paths.filename_script_graphviz

			logger.info("Generating the panels...")

			generator_plot_timeseries.plot(
				data_inference.table_genotypes,
				palette = current_palette.get_genotype_palette(),
				filename = filename_timeseries_genotypes
			)

			generator_panel_timeseries.run(
				timeseries_genotype = data_inference.table_genotypes,
				timeseries_trajectory = data_inference.table_trajectories,
				filename = filename_timeseries_panel,
				palette = current_palette,
			)
			logger.info("Generating the annotated muller plots...")
			generator_plot_muller.plot(
				muller_df = data_lineage.table_muller,
				color_palette = current_palette.get_genotype_palette(),
				annotations = genotype_annotations,
				filename = filename_mullerplot_annotated
			)

			generator_plot_muller_panel.plot(
				timeseries = data_inference.table_genotypes,
				muller_df = data_lineage.table_muller,
				palette = current_palette.get_genotype_palette(),
				annotations = genotype_annotations,
				filename = filename_muller_panel
			)

			logger.info("Generating the unannotated muller plots...")
			generator_plot_muller.plot(
				muller_df = data_lineage.table_muller,
				color_palette = current_palette.get_genotype_palette(),
				annotations = None,
				filename = filename_mullerplot_unannotated
			)
			logger.info("Generating the lineageplot...")
			lineageplot = graphics.flowchart(
				edges = data_lineage.table_edges.reset_index(),
				palette = current_palette.get_genotype_palette(),
				annotations = genotype_annotations,
				filename = filename_lineageplot,
				add_score = False
			)
			filename_graphviz_script.write_text(lineageplot.to_string())
			if suffix == 'png':
				# PIL doesn't work with svgs
				logger.info("Generating the lineage panel...")
				graphics.generate_lineage_panel(filename_mullerplot_annotated, filename_lineageplot,
					filename_lineage_panel)
