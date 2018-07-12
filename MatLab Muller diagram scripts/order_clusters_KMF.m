%SCRIPT: order_clusters
%DATE CREATED: 10-22-12
%
%CREATED BY: Katya Kosheleva
%
%PURPOSE: Given mutational clusters and their trajectories encoded in the 
%arrays 'allgenotypes' and 'gentrajectories', this program creates an array
%nesting these clusters, gennests(npop, gen, pnests), where npop is the 
%population index, gen is the index of the cluster, and pnest (in order) 
%are the clusters (i.e. mutations) occuring on the background of gen.
%
%GOAL: For each time point, give segregating genotypes, including
%backgrounds, (and another index denoting the probability of this 
%arrangement, including data from later times)  i.e.
%
%timepoint 1: {1, {1}, {2}}
%timepoint 2: {0.9, {1, 3}, {1, 4}, {2}}
%timepoint 3: {1, {1, 4, 5} }
%
%
% The way to do this is to sort the clusters by their maximum frequency 
% attained (since this guarantees that we go in order of nesting). Given a 
% cluster that we pick this way, there are two possibilities:
%(1) It occurs on the wild type background
%(2) It occurs on the background of an existing genotype
%
%Given (2), there are two possibilities: it nests farther into a
%genotype vector (i.e. {1, 4} -> {1, 4, 5}) or it splits an
%existing genotype vector in two ( {1, 4} -> {1, 4}, {1, 5}). 
%
%
%
%There are three criteria for nesting trajectories: 
%(1) If two trajectories, at a given time point, are at >50%, then either (1)
% one is on the background of the other, (2) they are actually part of the same
% cluster
%
%(2) Fixations: Everything that went extinct when something else fixed
% was in a different background
%
%(3) Derivatives: if increase of one is always (usually) associated with increase of
% another, then same background; if decrease, then different
% 
%(4) Low frequency genotypes associated according to probability (this is
%necessarily ambiguous)
%
%
%In the end, goal is to have a 14(+) column array, where (1) is the (full)
%genotype, (2) is the probability of this genotype, s.t. all genotypes
%containing every digit must sum to 1 at each time point, and 3-14 is the
%frequency of the genotype at each time point.
%

npops = 1;
frequencies = sort(0:.15:.90, 2 ,'descend');
deltacut = 0.01;
gensortedtotal = 0;
genneststotal = 0;

%EDITED BY: Kenneth Flynn
%added two new variables to allow for variable sized data sets based on
%user defined variables at run time.

trajSize = incolNum - 2;
nameSize = incolNum - 1;


for npop = 1
    trajectories = squeeze(gentrajectories(npop, any(squeeze(gentrajectories(npop, :, :)), 2),  :));  
    gensorted = [0; 0; 0];
    
    %%%%%%%%%%%%%STEP 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%SORT CLUSTERS BY MAXIMUM FREQUENCY ATTAINED%%%%%%%%%%%
    %%%%%%%IF MAXIMUM FREQUENCY ATTAINED IS 1, SORT BY TIME OF FIXATION%%%%
    
    %NOTE: FOR SORTING, FOR EVERY CLUSTER THAT GETS CAUGHT, LOOK +\- 10%
    %AND INCLUDE THOSE AS WELL
    
    %each row is a trajectory. for the ones that fixed or nearly fixed (> 85%), arrange by how
    %early they did so
    
    for time = 1:tNum
        if  ~isempty(trajectories(trajectories(:, time)>.85, : )) %if there are fixed trajectories at time
            indices = find(trajectories(:, time)>.85);
            newindices = indices(~ismember(indices, gensorted(1,:)));
            %             if size(newindices, 1) ==1 %if there's only one new fixing genotype at time
            %                 gensorted(1, size(gensorted, 2)+1) = newindices;
            %                 gensorted(2,  size(gensorted, 2)) = trajectories(newindices, time);
            %                 gensorted(3, size(gensorted, 2)) = time;
            %             else %otherwise arrange by size one time step before
            if ~isempty(newindices)
                temp = 0;
                for index = 1: size(newindices, 1) %arrange all fixed genotypes by when they were first detected and first above 15%
                    temp(index, 1) = newindices(index);
                    temp(index, 2) = find( trajectories(newindices(index), :) >0.03, 1, 'first');
                    temp(index, 3) = find( trajectories(newindices(index), :) >0.15, 1, 'first');
                end
                while size(temp,1)>1
                    [junk, where] = sort(temp(:, 2), 'ascend');
                    temp = temp(where,:); % sort according to values in second column from smallest to largest
                    if temp(1,2) - temp(2, 2) > 0 %pick out the one that was detected first
                        gensorted(1, size(gensorted, 2)+1) = temp(1, 1);
                        gensorted(2,  size(gensorted, 2)) = trajectories(temp(1, 1), time);
                        gensorted(3, size(gensorted, 2)) = time;
                        temp = temp(2:end,:); %remove this genotype from temp
                    else %otherwise pick the one that was > 15% first
                        [junk, where] = sort(temp(:, 3), 'ascend');
                        temp = temp(where,:);
                        gensorted(1, size(gensorted, 2)+1) = temp(1, 1);
                        gensorted(2,  size(gensorted, 2)) = trajectories(temp(1, 1), time);
                        gensorted(3, size(gensorted, 2)) = time;
                        temp = temp(2:end,:); %remove this genotype from temp
                    end
                end
                gensorted(1, size(gensorted, 2)+1) = temp(1, 1);
                gensorted(2, size(gensorted, 2)) = trajectories(temp(1, 1), time);
                gensorted(3, size(gensorted, 2)) = time;
            end
        end
    end
    
    if size(gensorted, 2)< size(trajectories, 1)+1 %if there are still unsorted trajectories, we basically do the same thing, but scale with frequency
        for frequency = 1:size(frequencies, 2)
            for time = 1:tNum
                if  ~isempty(trajectories(trajectories(:, time)>frequencies(frequency), : )) %if there are trajectories > frequency at time
                    indices = find(trajectories(:, time)>frequencies(frequency));
                    newindices = indices(~ismember(indices, gensorted(1,:)));
                    %             if size(newindices, 1) ==1 %if there's only one new fixing genotype at time
                    %                 gensorted(1, size(gensorted, 2)+1) = newindices;
                    %                 gensorted(2,  size(gensorted, 2)) = trajectories(newindices, time);
                    %                 gensorted(3, size(gensorted, 2)) = time;
                    %             else %otherwise arrange by size one time step before
                    if ~isempty(newindices)
                        temp = 0;
                        for index = 1: size(newindices, 1) %arrange all fixed genotypes by when they were first detected and first above 15%
                            temp(index, 1) = newindices(index);
                            temp(index, 2) = find(trajectories(newindices(index), :) >0.03, 1, 'first');
                            if ~isempty(find( trajectories(newindices(index), :) >0.15 , 1) )
                                temp(index, 3) = find( trajectories(newindices(index), :) >0.15, 1, 'first');
                            else
                                temp(index, 3) = 13;
                            end
                        end
                        while size(temp,1)>1
                            [junk, where] = sort(temp(:, 2), 'ascend');
                            temp = temp(where,:); % sort according to values in second column from smallest to largest
                            if temp(1,2) - temp(2, 2) > 0 %pick out the one that was detected first
                                gensorted(1, size(gensorted, 2)+1) = temp(1, 1);
                                gensorted(2,  size(gensorted, 2)) = trajectories(temp(1, 1), time);
                                gensorted(3, size(gensorted, 2)) = time;
                                temp = temp(2:end,:); %remove this genotype from temp
                            else %otherwise pick the one that was > 15% first
                                [junk, where] = sort(temp(:, 3), 'ascend'); 
                                temp = temp(where,:);
                                gensorted(1, size(gensorted, 2)+1) = temp(1, 1);
                                gensorted(2,  size(gensorted, 2)) = trajectories(temp(1, 1), time);
                                gensorted(3, size(gensorted, 2)) = time;
                                temp = temp(2:end,:); %remove this genotype from temp
                            end
                        end
                        gensorted(1, size(gensorted, 2)+1) = temp(1, 1);
                        gensorted(2,  size(gensorted, 2)) = trajectories(temp(1, 1), time);
                        gensorted(3, size(gensorted, 2)) = time;
                    end
                end
            end
        end
    end
    
    
     gensorted = gensorted(:, 2:end); %get rid of the initial 0.
     
     
     %%%%%%%%%%%%STEP 2%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%METHODICALLY GO THROUGH EACH CLUSTER IN THE LIST%%%%%%%
     %%%%%ACCORDING TO A LIST OF CRITERIA, ARRANGE IT IN THE CORRECT%%%%
     %%%%%BACKGROUND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
     gennests = zeros(1, incolNum);
     gennests(1, nameSize) = gensorted(1, 1); %name of first genotype
     gennests(1, 1) = 1; %probability that this is the real or correct genotype
     gennests(1, 2:trajSize) = trajectories(gensorted(1, 1),:); %trajectory of this genotype
     
     for type = 2:size(gensorted, 2)
         %COMPARE TO GENOTYPES GOING BACKWARDS IN ORDER, UNTIL 'YES'
         %HAPPENS
         %QUESTIONS:
                  
         delta = 0;
         
         for test = 1:(type-1)          
          % (1) is the sum of the frequencies consistently greater than 1?
         % if so, put it in the background of the other
             
             genotype = 0;
             testtrajectory = trajectories(gensorted(1, type-test), :);
             typetrajectory = trajectories(gensorted(1, type), :);
             trajectorysum = testtrajectory + typetrajectory;
            if size(trajectorysum(trajectorysum > 1.03), 2)>1 || ~isempty(trajectorysum(trajectorysum > 1.15)) %if the sum is larger than 1 at least twice, or much larger than 1 once
                    %SUCCESS!  
                     cut = gennests(type-test, nameSize: end);
                     genotype = [cut( cut ~= 0) gensorted(1, type) ] ;
                     gennests(type, 1) = 1; %probability that this is the real or correct genotype
                     gennests(type, 2:trajSize) = trajectories(gensorted(1, type),:); %trajectory of this genotype
                     gennests(type, nameSize: trajSize+ size(genotype, 2)) = genotype;
                     break
            end
                
         % (2) if the current genotype is ever >15% larger than the other,
         % then they are not on the same background
 
         trajectorydiff = testtrajectory- typetrajectory;
         if size(trajectorydiff(trajectorydiff < -0.02), 2)>1 || ~isempty(trajectorydiff(trajectorydiff < -0.15)) %if typetrajectory > testtrajectory more than twice, or much greater once             
             continue
         end
         
         %(3) Look at the first point when both trajectories are non-zero, and
         %observe Delta with time. Are derivatives correlated or
         %anti-correlated? If sufficiently correlated, they are on the same
         %background.
            
           startpoint = max( [find(testtrajectory > 0.02, 1, 'first' ) find(typetrajectory > 0.02, 1, 'first') ]);           
           endpoint = min( [find(testtrajectory > 0.02, 1, 'last' ) find(typetrajectory > 0.02, 1, 'last') ]); 
            
           deltatest = 0;
           deltatype = 0;
           
%           if endpoint-startpoint>2
               for time = startpoint:endpoint-1
                   %compare derivatives here
                   deltatest(time+1-startpoint) = testtrajectory(time+1) - testtrajectory(time);
                   deltatype(time+1-startpoint) =  typetrajectory(time+1) - typetrajectory(time);
               end
%               testsize = dot(deltatest, deltatest);
%               typesize = dot(deltatype, deltatype);
%               deltatest = deltatest/sqrt(testsize); %making these unit vectors
%               deltatype = deltatype/sqrt(typesize);
               delta(test) = dot(deltatype, deltatest);
               
               if delta(test) > deltacut %then they are (probably) on the same background
                   %however, need to do one last test- these two genotypes
                   %cannot sum to larger than the background
                   cut = gennests(type-test, nameSize: end);
                   genotype = [cut( cut ~= 0) gensorted(1, type) ] ;
                   gennests(type, 1) = 1; %probability that this is the real or correct genotype
                   gennests(type, 2:trajSize) = trajectories(gensorted(1, type),:); %trajectory of this genotype
                   gennests(type, nameSize: trajSize+ size(genotype, 2)) = genotype;
                   break
               end
               
               %EVENTUALLY: will need to use the value of delta to assign
               %probabilities for each background. Will have to keep
               %'significant' values of delta, and then weigh each
               %possibility by its likelihood
               
               %Finally, there is one more consideration: the sum of each
               %background at every time point cannot exceed 1               
               
%           end
         end
        
               if ~ismember(gensorted(1, type), gennests(:, nameSize:end)) %if it hasn't been matched with a background, two things can happen
                   %(1) IF it is logically possible (i.e., if the sum of it
                   %and all existing backgrounds at each time point is ~1 or less),
                   %then we can make it its own background. otherwise put
                   %it in with the genotype it's most correlated with
                   backgrounds = gennests( gennests(:,tNum) == 0 , 2:trajSize);
                   total = sum(backgrounds,1) + typetrajectory;
                    if size(total(total >1), 2)>1 || ~isempty(total(total < 1.15))
                        gennests(type, 1) = 1; %probability that this is the real or correct genotype
                        gennests(type, 2:trajSize) = trajectories(gensorted(1, type),:); %trajectory of this genotype
                        gennests(type, nameSize) = gensorted(1, type);
                    else
                        if ~isempty(delta(delta>0)) %if this cluster is positively correlated with something, stick it there
                            test = find(delta == max(delta));
                            cut = gennests(type-test, nameSize: end);
                            genotype = [cut( cut ~= 0) gensorted(1, type) ] ;
                            gennests(type, 1) = 0.5; %just a dummy value for now
                            gennests(type, 2:trajSize) = trajectories(gensorted(1, type),:); %trajectory of this genotype
                            gennests(type, nameSize: trajSize+ size(genotype, 2)) = genotype;
                        else
                            fprintf('SOMETHING HAS GONE HORRIBLY WRONG FOR CLUSTER %f IN POPULATION %f\n',gensorted(1, type),npop);
                        end
                    end                 
               end               
         
         gensortedtotal(npop, 1:size(gensorted, 1), 1:size(gensorted, 2)) = gensorted;
         genneststotal(npop, 1:size(gennests, 1), 1:size(gennests, 2)) = gennests;
         
         

         
         %these three conditions should take care of all genotypes that
         %reach a frequency of >50%, and their probabilities should all be
         %1.
        
         
         
     end
     
    
end

clear backgrounds cut delta deltacut deltatest deltatype endpoint frequencies...
    frequency gennests genotype gensorted index indices junk newindices ...
    startpoint temp test testtrajectory time total trajectories trajectorydiff...
    trajectorysum type typetrajectory where










