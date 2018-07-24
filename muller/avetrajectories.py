from muller.time_series_import import import_timeseries
from muller.get_genotypes import get_genotypes
from muller import variables
import pandas
def get_mean_genotypes(genotypes, timeseries):
	mean_genotypes = list()
	for genotype in genotypes:
		genotype_timeseries = timeseries[timeseries['Trajectory'].isin(genotype)]
		mean_genotype_timeseries = genotype_timeseries.mean()
		mean_genotype_timeseries.name = "|".join(map(str,genotype))
		mean_genotypes.append(mean_genotype_timeseries)

	mean_genotypes = pandas.DataFrame(mean_genotypes)
	mean_genotypes.pop('Trajectory')
	mean_genotypes.pop('Position')
	return mean_genotypes

if __name__ == "__main__":
	time_series, info = import_timeseries(variables.filename)
	genotypes = get_genotypes(time_series)

	mean_genotypes = get_mean_genotypes(genotypes, time_series)
