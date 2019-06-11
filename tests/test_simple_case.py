import pandas
import pytest

from muller.dataio import import_table
from muller.inheritance import order


@pytest.fixture
def lineage_workflow() -> order.LineageWorkflow:
	return order.LineageWorkflow(0.03, 0.97, 0.03, 0.03, 0.01)


@pytest.fixture
def three_genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	4	5
	genotype-C	0	0	0	0.3	0.7	1
	genotype-A	0	0	0	0.1	0.5	0.5
	genotype-B	0	0.1	0.15	0.03	0	0
	"""

	return import_table(string, index = 'Genotype')


@pytest.fixture
def five_genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	5	8	13	21	34
	genotype-C	0	0	0	0.3	0.7	1	1	1	1
	genotype-A	0	0	0	0	0.45	0.5	0.55	0.7	0.85
	genotype-E	0	0	0	0	0	0	0.05	0.55	0.5
	genotype-B	0	0.07	0.1	0.02	0.01	0	0	0	0
	genotype-D	0	0	0	0	0	0	0.07	0	0.01
	"""
	return import_table(string, index = 'Genotype')


@pytest.fixture
def ten_genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	5	8	13	21	34
	genotype-C	0.00	0.00	0.00	0.30	0.70	1.00	1.00	1.00	1.00
	genotype-J	0.00	0.00	0.00	0.00	0.70	0.80	1.00	1.00	1.00
	genotype-A	0.00	0.00	0.00	0.00	0.45	0.50	0.55	0.70	0.85
	genotype-E	0.00	0.00	0.00	0.00	0.00	0.00	0.05	0.55	0.50
	genotype-F	0.10	0.33	0.20	0.10	0.05	0.00	0.00	0.00	0.00
	genotype-B	0.00	0.07	0.10	0.02	0.01	0.00	0.00	0.00	0.00
	genotype-G	0.00	0.00	0.07	0.02	0.10	0.00	0.00	0.00	0.00
	genotype-H	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.10
	genotype-I	0.00	0.00	0.00	0.00	0.00	0.07	0.06	0.00	0.00
	genotype-D	0.00	0.00	0.00	0.00	0.00	0.00	0.07	0.00	0.01
	"""
	return import_table(string, index = 'Genotype')


def test_three_genotypes(three_genotypes, lineage_workflow):
	expected = {
		'genotype-A': 'genotype-C',
		'genotype-B': 'genotype-0',
		'genotype-C': 'genotype-0'
	}

	# clusters = order.order_clusters(three_genotypes, 0.03, 0.97, 0.03, 0.03, 0.01)
	clusters = lineage_workflow.run(three_genotypes)
	assert clusters.as_dict() == expected


def test_five_genotypes(five_genotypes, lineage_workflow):
	expected = {
		'genotype-A': 'genotype-C',
		'genotype-B': 'genotype-0',
		'genotype-C': 'genotype-0',
		'genotype-D': 'genotype-C',
		'genotype-E': 'genotype-A'
	}
	result = lineage_workflow.run(five_genotypes)
	# result = order.order_clusters(five_genotypes, 0.03, 0.97, 0.03, 0.03, 0.01)
	assert result.as_dict() == expected


def test_ten_genotypes(ten_genotypes, lineage_workflow):
	expected = {
		'genotype-A': 'genotype-J',
		'genotype-B': 'genotype-F',
		'genotype-C': 'genotype-0',
		'genotype-D': 'genotype-A',
		'genotype-E': 'genotype-A',
		'genotype-F': 'genotype-0',
		'genotype-G': 'genotype-0',
		'genotype-H': 'genotype-E',
		'genotype-I': 'genotype-A',
		'genotype-J': 'genotype-C'
	}
	# Manual override since the ancestry of these genotypes is somewhat ambiguous.
	expected['genotype-H'] = 'genotype-J'
	expected['genotype-I'] = 'genotype-C'
	expected['genotype-D'] = 'genotype-J'
	# result = order.order_clusters(ten_genotypes, 0.03, 0.97, 0.03, 0.03, 0.01)
	result = lineage_workflow.run(ten_genotypes)
	assert result.as_dict() == expected
