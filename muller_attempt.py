from muller.import_data import import_table_from_string

if __name__ == "__main__":
	string = """
		Trajectory	0	17	25	44	66	75	90
		9	0	0	0	0	0	0.269	0.34
		17	0	0	0	0	0	0.266	0.312
		19	0	0	0	0.188	0.171	0.232	0.244
	"""
	table = import_table_from_string(string, index = 'Trajectory')
	print(table.to_string())