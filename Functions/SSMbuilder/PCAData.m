function [ssmV,Eval,Evec,MEAN,PCcum]=PCAData(RR)

% PCA by SVD

%INPUT
%   RR : M x N matrix with M equal to 3 times the number of vertices and N the number
%              of training data
%
%OUTPUT
%   
% ssmV : the shape vectors, an M x N Matrix
%   Eval : Eigen values of the data
%   Evec : Eigen vectors of the data
%   MEAN : Mean shape vector
%   PCcum: cumulative variance explained by the SSM vectors

% Mostly inspired by previous work and coding by Dirk Jan Kroon



s=size(RR,2);
MEAN=sum(RR,2)/s;

% Define residuals
resiudals=(RR-repmat(MEAN,1,s))/ sqrt(s-1);

% SVD 
[U2,S2] = svds(resiudals,s);
Eval=diag(S2).^2;
Evec=bsxfun(@times,U2,sign(U2(1,:))); 


Eval2=sqrt(Eval);

for i=1:length(Eval)
    ssmV(:,i)=Eval2(i,1)*Evec(:,i);
   
end

PCcum=cumsum(Eval)./sum(Eval);
