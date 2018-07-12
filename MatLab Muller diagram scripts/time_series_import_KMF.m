% SCRIPT: time_series_import

%created by K. Kosheleva, 10/12/2012
%
%This script imports the excel file 'all_nuclear_mutations_with_gene_counts*'
%and converts it to a matlab variable
%
%Inputs: none (but make sure the file is in the correct folder)
%Outputs: Matrix 'timeseries' and 'info'
%The matrix 'timeseries' has the format (x,n) where n is the number of mutations
%found and x corresponds to columns with the following information:
%1- Population Number
%2- Trajectory Number
%3- Position on W303 Reference
%4:15- Value at various timepoints
%16 - First above 0
%17 - Next above 10
%18 - s
%19- Fixed?
%20; Tfix
%21- T above 0.05
%22- T below 0.05
%23- Big?
%24- Nonsynonymous
%25- Synonymous
%26- Gene Hit Count
%27- 2 or more
%28- 3 or more
%29 - Neutral?
%30- SGDID
% Finally, the info cell matrix contains string info:
%
%info(1,:) = population
%info(2,:) = chromosome
%info(3,:) = mutation type, e.g. snp, del, etc
%info(4,:)  = nucleotide substitution
%info(5,:) = name of gene hit
%
%%

%names = ls('all_nuclear_mutations_Bt2.xls')

rawcolNum = 12 + tNum;
incolNum = 3 + tNum;
sheet = sheets;

[data, garbage] = xlsread(names(1,:), sheets(1,:)); %importing excel data

[a,b] = size(data);

timeseries = zeros(incolNum, a-1);
info = cell(5, a-1);

timeseries(1:2,:) = data(2:end,1:2)'; % POPULATION AND TRAJECTORY
timeseries(3,:) = data(2:end,4)'; %W303 REFERENCE
timeseries(4:incolNum,:) = data(2:end,13:rawcolNum)'; %FREQUENCIES AT TIME POINTS... columns 13-15 starting at generation 0 in column N....19timeseries(4:10,:) = data(2:end,13:19)'
%timeseries(11:13,:) = data(2:end,20:22)'; %FIRST ABOVE/NEXT ABOVE/ S
%timeseries(14:18,:) = data(2:end, 26:30)'; %FIXED/TFIX/TABOVE/TBELOW/BIG
%timeseries(19:25,:) = data(2:end, 33:39)'; %ALL THE REST...

info(1,:) = garbage(2:end, 1);
info(2,:) = garbage(2:end, 4);
info(3,:) = garbage(2:end, 6);
info(4,:) = garbage(2:end, 7);
info(5,:) = garbage(2:end, 8);

%keep only important stuff
clear a b data garbage

