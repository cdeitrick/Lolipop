# Genotype Filters
Some of the detected trajectories or genotypes do not represent real mutations, and need to be removed from the analysis. Since all pre-existing variation is removed when a specific genotype sweeps and fixes, any mutation that is seen before, during, and after the timepoint where the sweep occured are likely due to measurement error.

Common situations where a trajectory or genotype should be removed:

1. A mutation appears before and after a genotype sweeps and removes all pre-existing variation. Since the mutation should have been removed when the genotype fixed
2. A mutation fixes immediately during the same timepoint it is detected (there are no intermediate values), then becomes undetected again.

An exception has been made for the case where a mutation appears both before and after a genotype sweeps, but is undetected at the timepoint where this occurs. This represents a mutation arising, being removed during a genotype sweep, then arising again.

Note that when a genotype is removed, all of the trajectories that comprise that genotype are removed and the clustering algorithm is re-run on the remaining dataset. 
