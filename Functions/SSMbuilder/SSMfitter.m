function [RMSerror,ReallignedV,transform,SSMfit,EstimatedModes]=SSMfitter(MEAN,Fmodel,ssmV,V,F,nrofshapemodes)

% This function fit a previous developped SSM to a given shape.

% INPUT

%   MEAN : Mean shape vector ; a M*1 vector
%   Fmodel : Triangulation Matrix of the model; a P by 3 matrix
%   ssmV : the shape vectors, an M x N Matrix
%   V : vertices of the given shape
%   F : triangulation matrix of the given shape
%   nrofshapemodes : the suggested nr of shape mode vectores to be used

% OUTPUT

%   RMSerror : rooted mean squared error after fitting of the SSM
%   ReallignedV : the new position of the shape defined by V following the
%               fitting procedure
%   transform : transformationmatrix applied to get from V to ReallignedV
%   SSMfit :    Predicted shape by fiiting of the SSM to V
%   EstimatedModes : Estimates shape modes of V to fit the model




%EXAMPLE

% load('EXAMPLE1.mat')
% [ssmV,Eval,Evec,MEAN,PCcum,Modes]=SSMbuilder(Xfemur,Yfemur,Zfemur);
% [RMSerror,ReallignedV,transform,SSMfit,EstimatedModes]=SSMfitter(MEAN,Fmodel,ssmV,V,F,nrofshapemodes)

[p,q]=size(ssmV);
s=p/3;

MEANX=MEAN(1:s);
MEANY=MEAN(s+1:2*s);
MEANZ=MEAN(2*s+1:3*s);

MeanParentbone=[MEANX,MEANY,MEANZ];


[error,ReallignedVtarget,transform]=rigidICP(MeanParentbone,V);


BTXX=ssmV(1:p/3,:);
BTXY=ssmV(p/3+1:2*(p/3),:);
BTXZ=ssmV(2*(p/3)+1:end,:);

errortemp(1,:)=0;
index=2;

estimate=MeanParentbone;
[errortemp(index,:),Reallignedtargettemp,estimate,EstimatedModes]=ICPmanu_allignSSM(ReallignedVtarget,MeanParentbone,estimate,BTXX,BTXY,BTXZ,nrofshapemodes);
display ('error')
d=errortemp(index,:)


while (abs(errortemp(index-1,:)-errortemp(index,:)))>0.000001
[errortemp(index+1,:),Reallignedtargettemp,estimate,EstimatedModes]=ICPmanu_allignSSM(Reallignedtargettemp,MeanParentbone,estimate,BTXX,BTXY,BTXZ,nrofshapemodes);
index=index+1;
d=errortemp(index,:)

end

error=errortemp(index,:);
SSMfit=estimate;
[d,ReallignedV,transform] = procrustes(Reallignedtargettemp,V,'scaling',0);

[projections]=project(SSMfit,Fmodel,ReallignedV,F);

 RMSerror=sqrt(mean((sqrt(sum((projections-SSMfit).^2,2)).^2)));



