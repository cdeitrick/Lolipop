
 The scripts currently use the two-step method by default, but this will be replaced by the hierarchical method after more testing is done.
 The hierarchical method and the desired distance metric can be selected using the `--method` and `--metric` options, respectively.

 # Two-Step clustering vs hierarchical clustering

 For the B1 testing population, the two-step method produces the following clusters:

![B1.twostep.trajectories](figures/B1.twostep.trajectories.edited2.png)

Each line is a single trajectory and is grouped with those of the same color. A major weakness of the two-step clustering method is its inability to accurately extract genotypes during timepoints when a large number of mutational trajectories are detected at similar frequencies. An example of this has been highlighted in the image above, where there are ten trajectories split among three genotypes.

Let's use the tan genotype (consisting of 6 trajectories) as an example. The trajectories labelled A, B, and C appear to be more correlated to each other than to the other three trajectories that have been grouped into this genotype. While the two-step method determined that there was insufficient evidence to group the trajectories into a separate genotype, a human observer may disagree.

Combining the dinomial distance metric with hierarchical clustering provides the following plot:
![B1.hierarchy.binomial.trajectories](figures/B1.hierarchy.binomial.trajectories.png)

This method splits the previous genotype up and groups A and B together (blue), while C (yellow) is paired with another trajectory which was not part of the original genotype, but shares a similar path as C.
