size = 25;
gentrajectories = 0;
rowLetter = 'A';

for iter = 1:size
    rowNum = num2str(iter);
    Range = [rowLetter rowNum];
    xlswrite(sprintf('%sgclusterOUT.xls', names(23:end-4)), allgenotypes(1,:,iter), 'Sheet1',Range);
end