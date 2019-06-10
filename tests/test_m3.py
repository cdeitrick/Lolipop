import pandas
import pytest
from muller import dataio
from muller.inheritance import order


@pytest.fixture
def genotypes() -> pandas.DataFrame:
	string = """
	Genotype	0	1	2	3	5	6	7	8	9	10
	genotype-5	0.021	0.145	0.264	0.587	0.9615	1.00	1.00	1.00	1.00	1.00
	genotype-14	0.00	0.00	0.173	1.00	0.946	1.00	1.00	1.00	1.00	1.00
	genotype-11	0.00	0.00	0.170	0.55	0.947	1.00	1.00	1.00	1.00	1.00
	genotype-9	0.427	0.241	0.081	0.113	0.00	0.00	0.00	0.00	0.00	0.00
	genotype-15	0.00	0.00	0.00	0.00	0.01	0.325	0.077	0.252	0.261	0.00
	genotype-16	0.00	0.00	0.00	0.00	0.00	0.24	0.067	0.22	0.221	0.33
	genotype-2	0.00	0.00	0.00	0.00	0.00	0.263	0.07	0.081	0.069	0.042
	genotype-12	0.00	0.02925	0.02425	0.16275	0.00	0.00	0.00	0.00	0.00	0.00
	genotype-10	0.00	0.00	0.009	0.007	0.013	0.04	0.027	0.00	0.00	0.00
	"""
	return dataio.import_table(string, index = 'Genotype')


def test_order(genotypes):
	expected = {
		'genotype-2':  'genotype-11',
		'genotype-5':  'genotype-0',
		'genotype-9':  'genotype-0',
		'genotype-10': 'genotype-11',
		'genotype-11': 'genotype-5',
		'genotype-12': 'genotype-14',
		'genotype-14': 'genotype-0',
		'genotype-15': 'genotype-11',
		'genotype-16': 'genotype-11'
	}
	lineage_workflow = order.LineageWorkflow(.03, .97, .03, .03, .01)
	result = lineage_workflow.run(genotypes)
	#result = order.order_clusters(genotypes, .03, .97, .03, .03, .01)

	assert result.as_dict() == expected
