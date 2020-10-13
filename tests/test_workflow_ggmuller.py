from pathlib import Path
from typing import *
from muller import dataio, workflows
import pandas
from loguru import logger
import pytest


def get_table_edges() -> pandas.Series:
	edges = """
		Identity	Parent
		genotype-7	genotype-0
		genotype-10	genotype-7
		genotype-11	genotype-10
		genotype-8	genotype-7
		genotype-9	genotype-8
		genotype-1	genotype-0
		genotype-13	genotype-11
		genotype-4	genotype-1
		genotype-2	genotype-11
		genotype-5	genotype-0
		genotype-12	genotype-10
		genotype-6	genotype-0
		genotype-3	genotype-0
	"""
	return dataio.import_table(edges).set_index('Identity')['Parent']


def get_table_genotypes() -> pandas.DataFrame:
	genotypes = """
		Genotype	0	30	45	60	90
		genotype-7	0	0.167	0.55	0.91	0.972
		genotype-10	0	0.02	0.265	0.9	0.97
		genotype-11	0	0	0.15	0.836	0.945
		genotype-8	0	0.11363636	0.45	0.735	0.86
		genotype-9	0	0.08	0.475	0.625	0.848
		genotype-1	0	0.3345	0.28125	0.095	0.045
		genotype-13	0	0	0	0	0.315
		genotype-4	0	0.005	0.24	0.03	0.0375
		genotype-2	0	0	0	0	0.245666666666667
		genotype-5	0	0.139	0	0	0
		genotype-12	0	0	0	0	0.11
		genotype-6	0	0.063	0	0	0.01
		genotype-3	0	0	0	0	0.01
	"""
	return dataio.import_table(genotypes)


def get_table_population_smoothed() -> pandas.DataFrame:
	populations_smoothed = """
		Identity	Generation	Population
		genotype-7	0	0
		genotype-7	30	5.336364
		genotype-7	45	10
		genotype-7	60	1
		genotype-7	90	1
		genotype-10	0	0
		genotype-10	30	0
		genotype-10	45	11.5
		genotype-10	60	6.40000000000001
		genotype-10	90	1
		genotype-11	0	0
		genotype-11	30	0
		genotype-11	45	15
		genotype-11	60	83.6
		genotype-11	90	63
		genotype-8	0	0
		genotype-8	30	3.363636
		genotype-8	45	1
		genotype-8	60	11
		genotype-8	90	1
		genotype-9	0	0
		genotype-9	30	8
		genotype-9	45	47.5
		genotype-9	60	62.5
		genotype-9	90	84.8
		genotype-1	0	0
		genotype-1	30	32.95
		genotype-1	45	4.125
		genotype-1	60	6.5
		genotype-1	90	1
		genotype-13	0	0
		genotype-13	30	0
		genotype-13	45	0
		genotype-13	60	0
		genotype-13	90	31.5
		genotype-4	0	0
		genotype-4	30	0.5
		genotype-4	45	24
		genotype-4	60	3
		genotype-4	90	3.75
		genotype-2	0	0
		genotype-2	30	0
		genotype-2	45	0
		genotype-2	60	0
		genotype-2	90	24.5666666666667
		genotype-5	0	0
		genotype-5	30	13.9
		genotype-5	45	0
		genotype-5	60	0
		genotype-5	90	0
		genotype-12	0	0
		genotype-12	30	0
		genotype-12	45	0
		genotype-12	60	0
		genotype-12	90	11
		genotype-6	0	0
		genotype-6	30	6.3
		genotype-6	45	0
		genotype-6	60	0
		genotype-6	90	1
		genotype-3	0	0
		genotype-3	30	0
		genotype-3	45	0
		genotype-3	60	0
		genotype-3	90	1
		genotype-0	0	100
		genotype-0	30	29.65
		genotype-0	45	0
		genotype-0	60	0
		genotype-0	90	0
	"""

	return dataio.import_table(populations_smoothed, keep_empty = True)


def get_table_population_unsmoothed(table_genotypes: pandas.DataFrame) -> pandas.DataFrame:
	table_populations_unsmoothed = table_genotypes.melt(id_vars = ['Genotype'], value_vars = [0, 30, 45, 60, 90])
	table_populations_unsmoothed.columns = ['Identity', 'Generation', 'Population']
	table_populations_unsmoothed = table_populations_unsmoothed.sort_values(by = ['Identity', 'Generation'])
	# Convert to percentages
	table_populations_unsmoothed['Population'] = table_populations_unsmoothed['Population'] * 100



	return table_populations_unsmoothed


def test_workflow_ggmuller():
	table_edges = get_table_edges()
	table_genotypes = get_table_genotypes()
	table_populations_smoothed = get_table_population_smoothed()
	table_populations_unsmoothed = get_table_population_unsmoothed(table_genotypes)
	table_genotypes = table_genotypes.set_index('Genotype')

	# Run the tests

	result_smoothed = workflows.run_workflow_ggmuller(table_genotypes, table_edges, 0.03, smooth_values = True)
	# To make the comparison a little simpler:
	table_smoothed = result_smoothed.table_populations

	assert table_smoothed['Identity'].tolist() == table_populations_smoothed['Identity'].tolist()
	assert table_smoothed['Generation'].tolist() == table_populations_smoothed['Generation'].tolist()
	# assert table_smoothed['Population'].tolist() == pytest.approx(table_populations_smoothed['Population'].tolist(), rel = 1E-3)
	assert table_smoothed['Population'].tolist() == pytest.approx(table_populations_smoothed['Population'].tolist())

	result_unsmoothed = workflows.run_workflow_ggmuller(table_genotypes, table_edges, 0.03, smooth_values = False)
	table_unsmoothed = result_unsmoothed.table_populations.sort_values(by = ['Identity', 'Generation'])

	# Drop Genotype-0 in order to make the comparison simpler, since it would not be added to the unsmoothed expected table.
	table_populations_unsmoothed = table_populations_unsmoothed[table_populations_unsmoothed['Identity'] != 'genotype-0']
	table_unsmoothed = table_unsmoothed[table_unsmoothed['Identity'] != 'genotype-0']

	assert table_unsmoothed['Identity'].tolist() == table_populations_unsmoothed['Identity'].tolist()
	assert table_unsmoothed['Generation'].tolist() == table_populations_unsmoothed['Generation'].tolist()
	# assert table_smoothed['Population'].tolist() == pytest.approx(table_populations_smoothed['Population'].tolist(), rel = 1E-3)
	assert table_unsmoothed['Population'].tolist() == table_populations_unsmoothed['Population'].tolist()