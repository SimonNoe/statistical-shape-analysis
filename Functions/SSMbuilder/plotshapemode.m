function plotshapemode(nr_of_mode,ssmV,MEAN,Fdata)

% function to plot principal components of shape variation at the +3 standard
% deviation, average and -3 SD. 

%INPUT
% nr_of_mode : desired nr of shape mode to plot
% ssmV : set of principal components of SSM
% MEAN: mean shape vector
% Fdata : triangulation matrix of training data 



[p,q]=size(ssmV);
s=p/3;
SD = 2;

MEANX=MEAN(1:s);
MEANY=MEAN(s+1:2*s);
MEANZ=MEAN(2*s+1:3*s);

MeanParentbone=[MEANX,MEANY,MEANZ];


clf
trisurf(Fdata,MeanParentbone(:,1),MeanParentbone(:,2),MeanParentbone(:,3),'Edgecolor','none');
hold
colormap bone


light
lighting phong;
set(gca, 'visible', 'off')
set(gcf,'Color',[1 1 0.88])
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
view(0,90)

%collumn of ssmV gives nth principal component

PC=MEAN+SD*ssmV(:,nr_of_mode);
PCplus(:,1)=PC(1:s,1)+2*(max(PC(1:s,1))-min((PC(1:s,1))));
PCplus(:,2)=PC(s+1:2*s,1);
PCplus(:,3)=PC(2*s+1:3*s,1);
trisurf(Fdata,PCplus(:,1),PCplus(:,2),PCplus(:,3),'Edgecolor','none');



PC=MEAN-SD*ssmV(:,nr_of_mode);
PCmin(:,1)=PC(1:s,1)-2*(max(PC(1:s,1))-min((PC(1:s,1))));
PCmin(:,2)=PC(s+1:2*s,1,1);
PCmin(:,3)=PC(2*s+1:3*s,1);
trisurf(Fdata,PCmin(:,1),PCmin(:,2),PCmin(:,3),'Edgecolor','none');
