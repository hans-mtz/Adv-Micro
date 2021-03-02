function [pi_hat,N_j,h_N] =RUM_63_KernelPihat(budgets,shares,income,med_income,X,budget_l)
%% Code Description: Create pi_hat
% This module estimates the probability that patch patch i is chosen from
% budget set j, which we call pi_hat.  
% Inputs: 
%   - budgets:      median-income normalized budget lines
%   - shares:       share of expenditure on each composite good.
%   - income:       consumer income proxied with expenditure.
%   - med_income:   polynomial of median income over years
%   - X:            Vector determining valid patches on each budget B_j
%   - budget_l:     budget length
% Outputs:
%   - pi_hat:    probability of choosing patch i from budget j
%   - N_j:       Number of households in each year j
%   - h_N:       bandwidth

%% Global parameters
% This function is passed through a parfor loop.  Therefore, we cannot
% include global variables.  Instead, I put the variables that we need as
% an argument.  

%% Pre-assign variables
N_j       = zeros(budget_l,1);      % number of people
N_patches = size(X,1);              % Number of patches
pi_hat    = zeros(N_patches,1);     % probability of choosing patch i from budget j.
h_N       = zeros(budget_l,1);      % Bandwidth
d{budget_l} = [];

%% Patch choice d and N.  
% We loop through each budget (year)
for jj = 1:budget_l
    % Number of people in year jj
    N_j(jj) = size(shares{jj},1);
    
    % Reset temporary variables
    patch = zeros(N_j(jj),budget_l);  % Patch that person n is on
    d{jj} = zeros(N_j(jj),1);             % Matches patch to X  
    
    % Implied quantity of each good demanded by person n:
    quantity = shares{jj}./kron(ones(N_j(jj),1),budgets(jj,:));
    % That is, q = shares / (p/e) = shares*e/p = (pq/e)*e/p.
    % where e is medium income in period jj and p is price in period jj.
    
    % For each person n, check whether their purchase decision is 
    % affordable or not in other budgets B_t.  This will imply a unique 
    % patch. 
    value = quantity*(budgets.');
    patch(value>1) = 1;
    patch(value<1) = -1;
    patch(:,jj) = 0;
    
    % Match patch with rows in matrix X, so we get the implied patch that
    % person n is on.
    [~, d{jj}] = ismember(patch,X,'rows'); 
end

%% Kernel Smooting pi_hat
for jj = 1:budget_l
    % Standard Dev. of sample
    hat_sigma  = std(income{jj});
    
    % Silverman's Rule-of-Thumb Bandwidth
    h_N(jj,1) = 1.06*hat_sigma*N_j(jj)^(-1/5);
    
    % Kernel Weight using Silverman's Rule-of-Thumb (See, for instance, Hansen pg342)
    K_weight = normpdf((income{jj} - med_income(jj,1))/h_N(jj,1),0,1);
    
    % Find all unique patches of X that people purchase from.  
    % We ignore patches not purchased from, since pi_hat is zero for those.
    uniq_patch = unique(d{jj}); 
    for ii = 1:size(uniq_patch)
        % Calculate pi_hat using formula in Yuichi's notes.  
        ind = d{jj} == uniq_patch(ii);
        pi_hat(uniq_patch(ii)) = sum(K_weight(ind),1)/sum(K_weight,1);
    end
end

end
 






