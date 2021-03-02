function [pi_hat] = RUM_11_Pihat(Choice, X)
%% Code Description
% This function takes in as inputs Choice indicating the choice by person i
% in budget j and the set of patches X.  The output if a frequency
% estimator for choices.
%
% Input:
%   Choice      I by J matrix of choices
%   X           I by J matrix of patches
%
% Output:
%   pi_hat      I by 1 vector of frequency of choices

% Determine number of patches and number of budgets
[I,J]   = size(X);
N       = size(Choice,1);

% Preset output
pi_hat = zeros(I,1);

% Get patches
patch_size = zeros(J,1);
patch_ind  = zeros(J,1);
for jj = 1:J
    patch_size(jj,1) =  size(find(X(:,jj) == 0),1);
    patch_ind(find(X(:,jj) == 0),1) = jj;
end
patch_size = [0;patch_size];

% For each patch, estimate frequency
for jj = 1:J
    cc = sum(patch_size(1:jj));
    for ii = 1:patch_size(jj+1)
        pi_hat( cc + ii) = size(find(Choice(:,jj) ==  ii),1)/N;
    end   
end

end
