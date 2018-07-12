%SCRIPT: ordered_cluster_plots
%DATE CREATED: 10-22-12
%
%CREATED BY: Katya Kosheleva
%
%PURPOSE: Given arrays of nested clusters and their trajectories, 
%this script plots each cluster and includes the genetic background of each
%cluster in the legend.
%
%INPUTS: Need array genneststotal (created by the order_clusters script), in
%the workspace
%

npops = 1;
ncolors = 30;
%timepoints = [0 102 150 264 396 450 540];

ColorSet = varycolor(ncolors);

for iter = 1
    
    close

    plots = tight_subplot(1,  1, [0.06, 0.02]);

    set(figure(1),'Position',[80 80 3000 2000])
    
 %   for npop = 8*(iter-1)+1: 8*iter
        
        gennests = squeeze( genneststotal(npop, any(squeeze(genneststotal(npop, :, :)), 2),:));
        
        
        axes(plots(npop - trajSize*(iter-1)))      
         colorstep = floor(ncolors/size(gennests,1));
        
         entry = cell(size(gennests, 1), 1);
         
        for traj = 1: size(gennests, 1)
                plots(npop - trajSize*(iter-1)) = plot( timepoints, gennests(traj, 2:trajSize), 'o-', 'Color', ColorSet(1+colorstep*(traj-1),:), 'LineWidth', 2); 
                temp = gennests(traj, nameSize:end);
                entry{traj} = num2str( temp(temp ~= 0));
                hold on            
        end
                
        
         l=legend(entry, 'Location', 'NorthWest');
         %legend hide
        
        if rem(npop, 4) ==1    
            ylabel('frequency');
        end
        
        
        if rem(npop, 4) ~= 1
            set(gca, 'YTickLabel', '')
        end
        
        %if rem(npop-1, 8) >= 4
            xlabel(sprintf('time in %s', xUnits));
        %end
        
        str = sprintf('Population %g ,  %s',npop,  char(info(1, find(timeseries(1,:) == npop, 1, 'first')))); %char(info(1, trajectories(1))));
        title(str)
        
        %The next few lines make it so that you can see subsequent
        %populations by clicking the mouse
        
        
      
 %   end
 
 print figure2.eps -depsc
 
        w = waitforbuttonpress;
        
        while 1
            if w==0
                break
            else continue
            end
            
        end
    close
    
    
   % saveas(figure(1), ['nestedgenotypes' num2str(iter) '.m'])
    
end





















