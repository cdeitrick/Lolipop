from setuptools import setup

setup(
	name = 'muller_diagrams',
	version = '0.4.1',
	packages = ['dataio', 'graphics', 'palettes', 'clustering', 'clustering.methods', 'clustering.metrics', 'inheritance', 'muller_output'],
	package_dir = {'': 'muller'},
	url = 'https://github.com/cdeitrick/muller_diagrams',
	license = 'MIT',
	author = 'chris deitrick',
	author_email = 'cld100@pitt.edu',
	description = 'A set of scripts to cluster mutational trajectories into genotypes and cluster genotypes by background',
	install_requires = ['pandas', 'loguru', 'scipy', 'matplotlib', 'pygraphviz', 'seaborn', 'numpy']
)
