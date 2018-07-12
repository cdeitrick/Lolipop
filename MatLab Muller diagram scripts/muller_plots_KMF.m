%SCRIPT: muller_plots
%
%CREATED BY: Katya Kosheleva
%
%PURPOSE: Given an array of frequency trajectories for each mutation/
%cluster of mutations, along with a nesting array for these mutations,
%this script plots the corresponding muller diagram.
%
%INPUTS:
%
%frequencies: m x t array of measured frequencies for each mutation/ 
%cluster of mutations.
%frequencies(i,j) gives the measured frequency of mutation cluster i at
%timepoint j.
%
%nests: m x m array of mutation clusters and their backgrounds. each row of
%the array gives the genetic background of a mutational cluster; i.e. a
%vector [1 3 5] denotes that mutation cluster 5 also carried mutation
%clusters 3 and 1 (note that the nesting array usually has trailing zeros
%after the last index, which do nothing).
%
%timepoints: vector of timepoints in generations corresponding to each
%trajectory in the frequency matrix.
%
%
%OUTPUTS: None, except for a plot of the muller diagram.
%

function[] = muller_plots(frequencies, nests, timepoints)

    transp = 0.7;

    m = size(frequencies, 1); %number mutations
    t = size(frequencies, 2); %number timepoints

    ncolors = 60;
    ColorSet = varycolor(ncolors);
    
    colorstep = floor(ncolors/size(nests, 1));

    maxnest = size(nests, 2);
    nests(:, size(nests,2)+1) = 0; %append a 0 to the end. just for convenience in loop plotting
    
    mullerplot = zeros(2, size(nests, 1), t);
    stack = zeros(t);
    space = zeros(t);
    previndices = 0;



    for nest = 1:maxnest
        indices = find(nests(:, nest+1) == 0 & nests(:, nest) ~= 0);
        if nest==1 %if there's no background to worry about
            for time = 1:t
                stack(time) = min(sum(frequencies(indices, time)), 1);
                space(time) = (1- stack(time))/(size(indices, 1)+1); %this is going to be the space between trajectories and the plot edges at a given timepoint
                totspace = 0; %running count of totspace makes sure that trajectories are layered on top of each other.
                scalefactor = min(1, stack(time)/sum(frequencies(indices, time))); %making sure all trajectories sum to 1
                for index = 1: size(indices, 1)
                    average = frequencies(indices(index), time)/2;
                    mullerplot(1, indices(index), time) = totspace+ space(time)+scalefactor*(average - frequencies(indices(index), time)/2);
                    mullerplot(2, indices(index), time) = totspace+  space(time)+scalefactor*(average + frequencies(indices(index), time)/2);
                    totspace = mullerplot(2, indices(index), time);
                end
            end
        else
            for previndex = 1:size(previndices, 1)
                subindices = find(nests(:, nest+1) == 0 & nests(:, nest) ~= 0 & nests(: , nest-1) == nests(previndices(previndex), nest-1));
                if ~isempty(subindices) %if one of these backgrounds has a trajectory inside it
                    for time = 1:t
                        prevalue = mullerplot(2, previndices(previndex), time)-  mullerplot(1, previndices(previndex), time);
                        stack(time) = min(sum(frequencies(subindices, time)), prevalue);
                        space(time) = (prevalue- stack(time))/(size(subindices, 1)+1);
                        totspace = 0;
                        scalefactor = min(1,  stack(time)/sum(frequencies(subindices, time))); %multiply everything by this scaling to make sure that subtrajectories
                        %'fit' inside their alleged background. this is to
                        %eliminate confusion caused by sampling uncertainty
                        %-> Perhaps should record this and note if any
                        %scalefactors are consistently too large?
                        for subindex = 1: size(subindices, 1)
                            average = frequencies(subindices(subindex), time)/2;
                            mullerplot(1, subindices(subindex), time) = totspace+ space(time)+mullerplot(1, previndices(previndex), time)+ scalefactor*(average - frequencies(subindices(subindex), time)/2);
                            mullerplot(2, subindices(subindex), time) = totspace+  space(time)+mullerplot(1, previndices(previndex), time)+scalefactor*(average + frequencies(subindices(subindex), time)/2);
                            totspace = mullerplot(2, subindices(subindex), time) - mullerplot(1, previndices(previndex), time) ;
                        end
                    end
                end
            end
        end

        previndices = indices; %record backgrounds in previous step

        hold on

        for index = 1:size(indices, 1)
            first = find(frequencies(indices(index), 1:t) > 0.001, 1, 'first'); %find first time it's nonzero
            last = find(frequencies(indices(index), 1:t) > 0.001, 1, 'last'); %find last time it's nonzero
            jbfill( timepoints(max(first-2, 1):min(t, last+2)), squeeze(mullerplot(2, indices(index), max(first-2, 1):min(t, last+2)))', squeeze(mullerplot(1, indices(index), max(first-2, 1):min(t, last+2)))','w', 'k', 1, 0.99 );
            jbfill( timepoints(max(first-2, 1):min(t, last+2)), squeeze(mullerplot(2, indices(index), max(first-2, 1):min(t, last+2)))', squeeze(mullerplot(1, indices(index), max(first-2, 1):min(t, last+2)))',ColorSet(1+colorstep*(indices(index)-1),:), 'k', 1, transp );
        end
        
    end
    
    axis([0 timepoints(end) 0 1])

    ylabel('frequency')
    xlabel('time in generations')
    set(gca,'Color',[ 1 1 1]);

    
    
end

