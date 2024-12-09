function [u0,up0] = Consistent_IC(J,M,state0,WS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% finds a consistant IC for the equation M*u'=J*u+b

% construct matrix
Mat1=[-J,M];
Mat2=speye(size(Mat1,2));
Mat=[Mat1;Mat2(1:size(Mat1,1),:)];
N=sum(WS.Nx);

% identify columns of all zeros
% this finds all boundary condition rows except 1
[~,col]=find(Mat);
col=unique(col);
inds=zeros(size(Mat,1)-length(col),1);
cnt=0;
for i=1:size(Mat,1)
    if ~any(i==col)
        cnt=cnt+1;
        inds(cnt)=i;
    end
end

for i=1:length(inds)
    Mat(inds(i),:)=0;
    Mat(inds(i),inds(i))=1;
end

b=[zeros(2*N,1);state0];
b(1)=-b(2*N+1);
b(2*N+1)=0;

phi=Mat\b;
u0=phi(1:2*N);
up0=phi(2*N+1:end);
end

