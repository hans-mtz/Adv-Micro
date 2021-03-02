function A=RUM_41_genA(X)
%% Code Description: Generate A
% This function generates A matrix from a given X matrix
% The columns of A represent pseudo-agents that have preferences consisten
% with GARP.  
% The rows of A  are 0 or 1 indicating the patch the pseudo-agent picks
% from.

% Globals
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores

%% Pre-define variables
[I,J]=size(X);
counter=sum(X==0);
A{counter(1)}= [];

%% Prebuild the revealed-preferred matrix 
% RP represents the directed graph of patches that are revealed preferred
% to other patches.  RP(i,j) = 1 iff patch i is revealed preferred to patch
% j.
RP = false(I,I);
row = 0;
for ii = 1:J
    RP(row + 1:row + counter(ii),:) = repmat((X(:,ii)<=0).',[counter(ii),1]);
    row = row+counter(ii);
end
 RP(eye(size(RP))~=false)=false;

%% Start tree-crawling algorithm
% Number of periods to check
A{counter(1)}= [];
parfor tt=1:counter(1)
    column=false(I,1);
    column(tt)=true;
    A{tt}=RUM_42_genAup(column,2,counter,RP,logical([]));
end

%% Reshape A
A = [A{1:counter(1)}];


end