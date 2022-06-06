function [X,Y,Z]=reallign2(X,Y,Z)

%iterative function realligns training data according to the mean and excluding
%outliers at the 0.05 level

[m,n] = size(X);
Template=horzcat(X(:,1),Y(:,1),Z(:,1));


for i=2:n
Ytemp= horzcat(X(:,i),Y(:,i),Z(:,i)); 
[outlierex]=outexallign(Template,Ytemp);
X(:,i)=outlierex(:,1);
Y(:,i)=outlierex(:,2);
Z(:,i)=outlierex(:,3);
end

Templatenew=horzcat(mean(X,2),mean(Y,2),mean(Z,2));
for i=1:n
Ytemp= horzcat(X(:,i),Y(:,i),Z(:,i)); 
[outlierex]=outexallign(Templatenew,Ytemp);
X(:,i)=outlierex(:,1);
Y(:,i)=outlierex(:,2);
Z(:,i)=outlierex(:,3);
end

Templatenew=horzcat(mean(X,2),mean(Y,2),mean(Z,2));
for i=1:n
Ytemp= horzcat(X(:,i),Y(:,i),Z(:,i)); 
[outlierex]=outexallign(Templatenew,Ytemp);
X(:,i)=outlierex(:,1);
Y(:,i)=outlierex(:,2);
Z(:,i)=outlierex(:,3);
end