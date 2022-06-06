function [error,Reallignedtarget,estimate,EstimatedModes]=ICPmanu_allignSSM(vnew,MEAN3d,estimate,BTXX,BTXY,BTXZ,nrofshapemodes)

[IDX1(:,1),IDX1(:,2)]=knnsearch(vnew,estimate);
[IDX2(:,1),IDX2(:,2)]=knnsearch(estimate,vnew);

IDX1(:,3)=1:length(estimate(:,1));
IDX2(:,3)=1:length(vnew(:,1));

Datasetsource=vertcat(MEAN3d(IDX1(:,3),:),MEAN3d(IDX2(:,1),:));
Datasettarget=vertcat(vnew(IDX1(:,1),:),vnew);

MEANaugmented=reshape(Datasetsource,size(Datasetsource,1)*3,1);

BTXaugmented=vertcat(BTXX(IDX1(:,3),:),BTXX(IDX2(:,1),:),BTXY(IDX1(:,3),:),BTXY(IDX2(:,1),:),BTXZ(IDX1(:,3),:),BTXZ(IDX2(:,1),:));


Zssm=reshape(Datasettarget-Datasetsource,size(BTXaugmented,1),1);
PI=pinv(BTXaugmented(:,1:nrofshapemodes));
EstimatedModes=PI*Zssm;
BTX=[BTXX(:,1:nrofshapemodes);BTXY(:,1:nrofshapemodes);BTXZ(:,1:nrofshapemodes)];
tempestimate=reshape(MEAN3d,size(MEAN3d,1)*3,1)+BTX*EstimatedModes;
estimate=reshape(tempestimate,size(tempestimate,1)/3,3);

IDX1=[];
IDX2=[];
[IDX1(:,1),IDX1(:,2)]=knnsearch(vnew,estimate);
[IDX2(:,1),IDX2(:,2)]=knnsearch(estimate,vnew);

IDX1(:,3)=1:length(estimate(:,1));
IDX2(:,3)=1:length(vnew(:,1));


m1=mean(IDX2(:,2));
s1=std(IDX2(:,2));
IDX1=IDX1(IDX1(:,2)<(m1+1.96*s1),:);

Datasetsource=vertcat(estimate(IDX1(:,3),:),estimate(IDX2(:,1),:));
Datasettarget=vertcat(vnew(IDX1(:,1),:),vnew);

[error,Reallignedtarget1] = procrustes(Datasetsource,Datasettarget,'scaling',0);
Reallignedtarget=Reallignedtarget1(size(IDX1,1)+1:end,:);
