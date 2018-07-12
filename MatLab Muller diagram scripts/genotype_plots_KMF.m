%SCRIPT: genotype_plots
%DATE CREATED: 10-19-12
%CREATED BY: Katya Kosheleva
%

%top: all trajectories
%middle: clustered trajectories
%bottom: average trajectory of each cluster


ncolors = 100;
ColorSet = varycolor(ncolors);

set(figure,'Position',[80 80 700 2000])


for npop = 1:npops
    
    clf
     
    plots = tight_subplot(3,  1,  0.02); 
    
    axes(plots(1))
    
    trajectories = timeseries(4:incolNum, timeseries(1,:) == npop);
    ntrajectories = size(trajectories, 2); %number trajectories
    
    colorstep1 = floor(ncolors/ntrajectories);
    p=0;
    
    hold off
    for i  = 1:ntrajectories
        p(i) = plot(timepoints, trajectories(:, i), '-o', 'Color', ColorSet(colorstep1*(i-1)+1,:));
        hold on
    end
    
    set(gca,'XTickLabel',[])
    
    title(int2str(npop))
    
    ylabel('frequency')
    %l = legend(p);
    
    genotypes = squeeze(allgenotypes(npop,any(squeeze(allgenotypes(npop, :, :)), 2),any(squeeze(allgenotypes(npop, :, :)), 1) ));
    
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
    
    set(gca,'XTickLabel',[])
    
    ylabel('frequency')
    
    %l = legend(p, 'Location', 'NorthWest');
    
    
    p=0;
    axes(plots(3))
    
    for i  = 1:ngenotypes
        p(i) = plot(timepoints, squeeze(gentrajectories(npop, i, :)), '-o', 'Color', ColorSet(colorstep2*(i-1)+1,:));
        hold on
    end
    
    ylabel('frequency')
    
    xlabel(sprintf('time (%s)', xUnits))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%END PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%PRINT%%%
    
   print figure1.eps -depsc
    
   w = waitforbuttonpress;
    
    while 1
        if w==0
            break
        else continue
        end
        
    end
    
    close


end

clear ColorSet colorstep1 colorstep2 genotypes i ncolors ngenotypes ...
    ntrajectories p plots t traj trajectories w