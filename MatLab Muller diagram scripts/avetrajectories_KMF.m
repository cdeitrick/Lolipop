%SCRIPT: avetrajectories
%DATE CREATED: 10-19-11
%CREATOR: KATYA KOSHELEVA
%
%quick script to extract average trajectories of a given genotype.
%Does not account for different read depths of mutations. all
%these variables are created/defined in the scripts time_series_import and
%get_genotypes. 
%
%gentrajectories is the output
%
%gentrajectories gives weighted average of each trajectory at each
%timepoint. 
%

npops = 1;
gentrajectories = 0;

for npop = 1: npops
    genotypes = squeeze(allgenotypes(npop,any(squeeze(allgenotypes(npop, :, :)), 2),any(squeeze(allgenotypes(npop, :, :)), 1) ));
    gens = size(genotypes, 1);
    for gen =  1:gens
            trajs = size(genotypes(gen, genotypes(gen,:) ~= 0), 2);
            trajectories = zeros(tNum, trajs); %changed 7 to tNum
            for traj = 1:trajs
                trajectories(:, traj) = timeseries(4:incolNum, timeseries(1,:) == npop & timeseries(2,:) == allgenotypes(npop, gen, traj));
            end
            for time = 1:tNum %changed 7 to tNum
                     gentrajectories(npop, gen, time) = sum(trajectories(time, :))/trajs;
            end
    end
end

%clear unnecessary variables
clear traj trajs trajectories gens gen genotypes

