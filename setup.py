from setuptools import setup
from muller import commandline_parser
setup(
	name = 'muller',
	version = commandline_parser.__VERSION__,
	packages = [
		'muller', 'muller.clustering', 'muller.dataio', 'muller.graphics', 'muller.inheritance',
		'muller.muller_output', 'muller.palettes', 'muller.clustering.metrics', 'muller.clustering.methods'
	],
	provides = 'muller',
	url = 'https://github.com/cdeitrick/muller_diagrams',
	license = 'MIT',
	author = 'chris deitrick',
	author_email = 'chrisdeitrick1@gmail.com',
	description = 'A set of scripts to cluster mutational trajectories into genotypes and cluster genotypes by background',
	install_requires = ['pandas', 'loguru', 'scipy', 'matplotlib','graphviz', 'pygraphviz', 'seaborn', 'numpy', 'xlrd'],
	tests_requires = ['pytest'],
	classifiers = [
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	scripts = ["muller/muller_workflow.py"]
)
