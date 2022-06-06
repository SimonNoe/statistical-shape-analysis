function [outlierex]=outexallign(doel,tedoen)

% exclusion of outliers

[d,tedoen,transform]=procrustes(doel,tedoen,'Scaling',0);

for i=1:20;
d=(doel-tedoen);
distance = sqrt(sum(d.^2, 2));
M=mean(distance);
ss=std(distance);
distance(:,2)=1:length(distance);
distance=distance(distance(:,1)<M+1.65*ss,:);
testdoel=doel(distance(:,2),1:3);
testtedoen=tedoen(distance(:,2),1:3);
[dnieuw,localnieuw,transform]=procrustes(testdoel,testtedoen,'Scaling',0);
tedoen=tedoen*transform.T + repmat(transform.c(1,1:3),length(tedoen(:,1)),1);
end

outlierex=tedoen;


