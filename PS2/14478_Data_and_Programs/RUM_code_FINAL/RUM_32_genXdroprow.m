function [ drop_row ] = RUM_32_genXdroprow(budgets,X,dimx,ii,J)
%% Code Description: RUM_33_genXdroprow
% This code determines whether or not we drop a particular row of X.  Due
% to problems with parfor, genX module needs to call this one rather than
% embed it inside genX module.

warning('off','MATLAB:nargchk:deprecated')

A = [];
b = [];
for tt=1:J %%translate (0,-1,1) into linear inequalities
    if X(ii,tt)==0
        A=[A;budgets(tt,:);-budgets(tt,:)];
        b=[b;1;-1];
    elseif X(ii,tt)==1
        A=[A;-budgets(tt,:)];
        b=[b;-1];
    else
        A=[A;budgets(tt,:)];
        b=[b;1];
    end
end
A=[A;-eye(dimx)];
b=[b;zeros(dimx,1)];
l=length(b);
% Pass A,b,l,dimx onto CVX.  If the objective function is 0, then there
% is no solution to Ax <= b. We drop that patch of X.
cvx_begin quiet;
   variable exes(dimx);
   bound=zeros(l,1);
   A*exes-b<=bound; 
cvx_end;
if cvx_optval==+Inf
    drop_row = 1;
else
    drop_row = 0;
end
end
