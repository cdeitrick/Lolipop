%SCRIPT: get_genotypes
%created by K. Kosheleva, 10/19/2012
%
%This script takes the matlab variable 'timeseries', and sorts all
%trajectories based on genotypes.
%Results corresponding to one genotype are sent to a new array
%containing the alleles in the genotype and their frequencies.
%
%Inputs: 
% (1) make sure that 'timeseries' is in the workspace
% (2) npops = number of populations in the dataset
% (3) linkcut = minimum probability relation between two trajectories for them to be related
% (4) relcut = any two trajectories in a genotype must be related by at least this much, or the genotype gets split in two
% (5) ncolors = number of colors to vary
%
%Outputs: 'allgenotypes': array of form (npop, genotype, trajectory) sorting trajectories in
%population npop according to the index of their genotype
%
%For now I ignore the different read depths of called mutations. 
%
%
%Edited by KBH 05/15/18
%matlab version 2018a does not have the combntns function. I replaced this function with nchoosek, which is what matlab suggests for compatability. Checked with 3 data sets and all appear to do the same analysis as the previous version of this script.


npops = 1;
plotting = 0;

%timepoints = [0 102 150 264 396 450 540]; % defined outside

relcut = 0.1; 
linkcut = 0.05; 
ncolors = 60;
totals = 0;

ColorSet = varycolor(ncolors);

allgenotypes = zeros(npops, 25, 25);
parray = 0;

for npop = 1
    trajectories = timeseries(4:incolNum, timeseries(1,:) == npop); %changed 10 to incolNum
    ntrajectories = size(trajectories, 2); %number trajectories
    combos = nchoosek(1:ntrajectories,2);
    for pair = 1:nchoosek(ntrajectories, 2)
        test = trajectories(:, combos(pair,:));
        test = test(test(:, 1)<0.97 | test(:, 2)<0.97, :);  %cut down to those where at least one is not fixed
        test = test(test(:, 1)>0.03 | test(:, 2)>0.03, :);
        ps = (test(:,1) + test(:, 2))/2; %average frequency between the two points
        sigma2 = prod([ps'; (1-ps)'])/5;
        dif = test(:, 1) - test(:,2); %these are consistent with n independent draws from a normal distribution,
        %assuming an average variance of  n_{binom} p (1-p), where p is the
        %average of the two frequencies, and n_{binom} is picked arbitrarily as a
        %value that gives reasonable uncertainties given our data
        %
        %n random draws will have
        %\sigma_{tot}^2 = \sum_i \sigma_i^2 =  n_{binom} \sum_i p_i(1 - p_i), from
        %a property of normal distributions
        %
        %for \bar{X} = \sum_{i}^{n_X} X/n_X
        %
        %\sigma_{\bar{X}} = \sigma_{tot}/n_X = \sqrt( \sum_i p_i(1-p_i))/\sqrt{n_X}
        %
        %Finally, given \sigma and \bar{X} we can construct a p-value for
        %the measurement by numerically computing the following integral:
        %
        %1 - \int_{- \bar{X}/ \sigma_{\bar{X}}}^{\bar{X}/
        %\sigma_{\bar{X}}} \frac{1}{\sqrt{2 \pi}} e^{-x^2/2} dx
        %
        
        
        
        
        sigmapair = sqrt(sum(sigma2))/size( dif, 1);
        difbar = sum(abs(dif))/size( dif, 1);
        
        pval = 1-erf(difbar/(sqrt(2)*sigmapair)); %probability that we are in the same genotype
           
        
        parray(npop, pair, 1:2) = combos(pair, :);
        parray(npop, pair, 3) = pval;
        parray(npop, pair, 4) = sigmapair;
        parray(npop, pair, 5) = difbar;
        
    end
    
    genotypes = 1; %by default the first trajectory makes a genotype category. rows are genotypes, columns members
    
    for pair = 1:nchoosek(ntrajectories, 2)
        if parray(npop, pair, 3) > relcut %are the genotypes related?
            if ismember(combos(pair, 1), genotypes) %is one of the trajectories (1) already listed?
                [row1, column1] = find(genotypes == combos(pair, 1));
                if ismember(combos(pair, 2), genotypes) %is (2) listed?
                    [row2, column2] = find(genotypes == combos(pair, 2));
                    if row2 ~= row1 %if (1)/ (2) listed under different genotypes, combine the two. Otherwise do nothing
                        n1 = size(genotypes(genotypes(row1,:) ~= 0), 2); %number of trajectories in first genotype
                        n2 = size(genotypes(genotypes(row2, :) ~= 0), 2); %number of trajectories in second genotype
                        genotypes(row1,(n1+1):(n1+n2)) = genotypes(row2, 1:n2); %add all members of genotype 2 to end of genotype 1
                        genotypes = genotypes(~ismember( 1:size(genotypes, 1), row2),:); %get rid of genotype 2
                    end
                else %if (1) listed and not (2), add (2) to genotype (1)
                    n1 = size(genotypes(genotypes(row1,:) ~= 0), 2);
                    genotypes(row1,n1+1) = combos(pair, 2);
                    %if neither are listed, put both in a new genotype with
                    %the two of them together
                end
            else %(1) not listed
                if  ismember(combos(pair, 2), genotypes) %(2) listed -> add (1) to genotype (2)
                    [row2, column2] = find(genotypes == combos(pair, 2));
                    n2 = size(genotypes(genotypes(row2,:) ~= 0), 2);
                    genotypes(row2,n2+1) = combos(pair, 1);
                else %neither listed
                    newrow = size(genotypes, 1);
                    genotypes(newrow+1, 1) = combos(pair, 1);
                    genotypes(newrow+1, 2) = combos(pair, 2);
                end
            end
        end
    end
    
    %at the end, look at all trajectories that are not listed and
    %append them as their own category.
    for traj = 1:ntrajectories
        if ismember(traj, genotypes) == 0
            newsize = size(genotypes, 1);
            genotypes(newsize+1, 1) = traj;
        end
    end
    
    %finally, for each genotype, make sure each trajectory pair has some
    %non-trivial linkage (say, >0.0005). this avoids falsely linking together
    %trajectories. if not, divide the offending trajectories into two
    %camps, and sort the rest according to which one they are more
    %closely linked to. repeat until everything is linked.
    
      ngenotypes = size(genotypes, 1)-1;
      
      while ngenotypes<size(genotypes,1) %while each loop continues to split genotypes
          
          ngenotypes = size(genotypes, 1);
          
          for geno = 1:ngenotypes
              pgens = 0;
              gens = genotypes(geno, genotypes(geno,:) ~= 0);
              if size(gens, 2) > 1
                  gencombos = nchoosek(gens,2);
                  for c = 1:size(gencombos,1) %getting pvalues between each combo
                      gencombos(c, 3) = parray(npop, ismember(squeeze(parray(npop,:,1:2)), [gencombos(c,1) gencombos(c,2)], 'rows') |  ismember(squeeze(parray(npop,:,1:2)), [gencombos(c,2) gencombos(c,1)], 'rows'), 3);
                  end
                  if all(gencombos(:,3)>linkcut)==0 %if some of the trajectories are unlinked
                      %badpairs = gencombos(gencombos(:,3) <= linkcut,:); %these are all the unlinked trajectories.
                      
                      [t, t] = min(abs(gencombos(:, 3) - linkcut));
                      %try picking two more closely linked trajectories
                      
                      newtype1 = gencombos(t,1); %assign one of the unlinked to a new genotype
                      newtype2 = gencombos(t, 2); %and the other unlinked to the other one
                      for gen = 1:size(gens, 2)
                          debug_element = gens(gen);
                          if ~ismember(gens(gen), [newtype1, newtype2]) %if the trajectory is not yet sorted
                              p1 = gencombos(ismember(gencombos(:, 1:2), [newtype1(1,1) gens(gen)] ,'rows') | ismember(gencombos(:, 1:2), [gens(gen) newtype1(1,1)] ,'rows'),3);
                              p2 = gencombos(ismember(gencombos(:, 1:2), [newtype2(1,1) gens(gen)] ,'rows') | ismember(gencombos(:, 1:2), [gens(gen) newtype2(1,1)] ,'rows'),3);
                              if p1 >= p2
                                  newtype1(1, size(newtype1, 2)+1) = gens(gen);
                              else
                                  newtype2(1, size(newtype2, 2)+1) = gens(gen);
                              end
                          end
                      end
                      %now that all of our old genotypes are reorganized, replace
                      %the old genotype
                      genotypes(geno,:) = 0;
                      genotypes(size(genotypes,1)+1, 1:size(newtype1, 2)) = newtype1;
                      genotypes(size(genotypes, 1)+1, 1:size(newtype2, 2)) = newtype2;
                  end
              end
          end
      end
 
          
          %finally, cut down the genotype array by removing rows with all
          %zeros:
           genotypes = genotypes(any(genotypes, 2), :);
           totals = totals + size(genotypes,1);
          
           allgenotypes(npop, 1: size(genotypes, 1), 1:size(genotypes, 2)) = genotypes;
 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plotting
    
        hold off

        clf
        %figure

        plots = tight_subplot(2,  1,  0.1);
        axes(plots(1))

        colorstep1 = floor(ncolors/ntrajectories);
        p=0;
        for i  = 1:ntrajectories
            p(i) = plot(timepoints, trajectories(:, i), '-o', 'Color', ColorSet(colorstep1*(i-1)+1,:));
            hold on
        end

        title('all trajectories')
        
        l = legend(p);

        ngenotypes = size(genotypes, 1);

        colorstep2 = floor(ncolors/ngenotypes);
        p=0;
        axes(plots(2))

        for i  = 1:ngenotypes
            traj = genotypes(i, genotypes(i,:) ~= 0);
            for t = 1:size(traj, 2)
                p(i) = plot(timepoints, trajectories(:, traj(t)), '-o', 'Color', ColorSet(colorstep2*(i-1)+1,:));
                hold on
            end
        end

        title('clustered trajectories')
        
        l = legend(p, 'Location', 'NorthWest');



        w = waitforbuttonpress;
        
        while 1
            if w==0
                break
            else continue
            end

        end
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%END PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
end
 

%some manual corrections to all genotypes:
%gen_corrections

%extract average trajectories from counts
%FOR DEBUGGING PURPOSES
avetrajectories_KMF

clear extraneous variables.

clear ans c column1 column2 combos dif difbar gen gencombos geno genotypes...
    gens n1 n2 ncolors newrow newsize newtype1 newtype2 ngenotypes npop ntrajectories...
    p1 p2 pair parray pgens plotting ps pval row1 row2 sigma2 sigmapair...
    t test time totals traj trajectories ColorSet linkcut relcut
