from setuptools import setup
from muller import commandline_parser
from pathlib import Path
if commandline_parser.DEBUG:
	message = "The scripts are currently in debug mode!"
	raise ValueError(message)
FOLDER = Path(__file__).parent
README = FOLDER / "README.md"

with README.open() as readmefile:
    LONG_DESCRIPTION = readmefile.read()

setup(
	name = 'muller',
	version = commandline_parser.__VERSION__,
	packages = [
		'muller', 'muller.clustering', 'muller.dataio', 'muller.graphics', 'muller.inheritance',
		'muller.muller_output', 'muller.palettes', 'muller.clustering.metrics', 'muller.clustering.methods'
	],
	extras_require = {
		'Show progressbar for large datasets': ['tqdm'],
		'Additional support for parsing files': ['beautifulsoup4']
	},
	provides = 'muller',
	url = 'https://github.com/cdeitrick/muller_diagrams',
	license = 'MIT',
	author = 'chris deitrick',
	author_email = 'chrisdeitrick1@gmail.com',
	description = 'A set of scripts to cluster mutational trajectories into genotypes and cluster genotypes by background',
	long_description = LONG_DESCRIPTION,
	long_description_content_type='text/markdown',
	install_requires = ['pandas', 'loguru', 'scipy', 'matplotlib','graphviz', 'pygraphviz', 'seaborn', 'numpy', 'xlrd'],
	tests_requires = ['pytest'],
	classifiers = [
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	scripts = ["Muller.py"]
)
