import pytest
import pandas
from import_data import import_table_from_string
from inheritance import order
@pytest.fixture
def sorted_df()->pandas.DataFrame:
	string = """
		Genotype	0	17	25	44	66	75	90
		genotype-1   0.0  0.00  0.261  1.0000  1.0000  1.0000  1.0000
		genotype-6   0.0  0.00  0.000  0.2730  0.7810  1.0000  1.0000
		genotype-3   0.0  0.00  0.000  0.0000  0.0000  1.0000  1.0000
		genotype-4   0.0  0.00  0.000  0.5250  0.4540  0.9110  0.9100
		genotype-5   0.0  0.00  0.000  0.1470  0.4500  0.9240  0.8870
		genotype-11  0.0  0.00  0.000  0.0000  0.2780  0.8220  0.8030
		genotype-2   0.0  0.38  0.432  0.0000  0.0000  0.0000  0.0000
		genotype-8   0.0  0.00  0.000  0.4030  0.4890  0.0570  0.0800
		genotype-14  0.0  0.00  0.000  0.0000  0.0000  0.2675  0.3260
		genotype-10  0.0  0.00  0.000  0.1380  0.2950  0.0000  0.0810
		genotype-12  0.0  0.00  0.000  0.0000  0.2335  0.1330  0.0375
		genotype-7   0.0  0.00  0.000  0.1880  0.1710  0.2320  0.2440
		genotype-9   0.0  0.00  0.117  0.0000  0.0000  0.0000  0.1030
		genotype-13  0.0  0.00  0.033  0.1060  0.1065  0.0000  0.0000
		genotype-15  0.0  0.00  0.000  0.1145  0.0000  0.1205  0.0615
	"""
	return import_table_from_string(string, index = 'Genotype')

def test_inheritance(sorted_df):
	pass