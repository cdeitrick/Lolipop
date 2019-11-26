from pathlib import Path

if __name__ == "__main__":
	filename = Path(__file__).parent / "tests" / "data" / "tables_ggmuller" / "B1_muller_try1.ggmuller.populations.tsv"
	import pandas
	table = pandas.read_csv(filename, sep = '\t')

	print(table.head().to_string())
	print("\n"*3)

	result = table.pivot(index = 'Identity', columns = 'Generation', values = 'Population')
	print(result.head().to_string())