function [ssmV,Eval,Evec,MEAN,PCcum,Modes]=SSMbuilder(Xdata,Ydata,Zdata)

%Framework to develop shape models. Allignment is performed with exclusion
%for outliers at the 0.05 level

%principal components of shape variation can be plotted using "plotshapemode.m"

% INPUT

% Xdata, Ydata en Zdata represent the x,y,z coordinates of the trainingdataset.
%               Each matrix is a M x N matrix with M equal to the number of vertices and N the number
%              of training data. 
%                 Data should be obtained by nonrigid registration of
%              a uniform template mesh. (see also nonrigidICP for registration and
%              patch remesher for template generation)

% OUTPUT

%   ssmV : the shape vectors, an M x N Matrix
%   Eval : Eigen values of the data
%   Evec : Eigen vectors of the data
%   MEAN : Mean shape vector
%   PCcum: cumulative variance explained by the SSM vectors
%   Modes: The shapemodes of the trainingdata within the shape model


%EXAMPLE

% load('EXAMPLE1.mat')
% [ssmV,Eval,Evec,MEAN,PCcum,Modes]=SSMbuilder(Xfemur,Yfemur,Zfemur);
%  plotshapemode(2,ssmV,MEAN,Ffem)
% With only 7 trainingdata and no optimisation of the registration using
% prior knowledge of the developped SSM, the results of present example can
% not be used for clinical or any other studies...

% allignment  of data
[Xdata,Ydata,Zdata]=reallign2(Xdata,Ydata,Zdata);

% construct of trainingdata matrix
RR=[Xdata;Ydata;Zdata];

% SVD decomposition (with thanks to previous work by Dirk Jan Kroon)
[ssmV,Eval,Evec,MEAN,PCcum]=PCAData(RR);

%definition of modes
[p,q]=size(RR);

ssmV(:,q)=[];

PI=pinv(ssmV);

for i = 1:q    
Scaling=[Xdata(:,i);Ydata(:,i);Zdata(:,i)]-MEAN;
Modes(i,:)=PI*Scaling;
end




