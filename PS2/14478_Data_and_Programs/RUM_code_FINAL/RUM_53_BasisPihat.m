function [pi_hat,N] =RUM_53_BasisPihat(budgets,shares,income,med_income,X,budget_l,poly_degree)
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
%   - poly_degree:  polynomial degree
% Outputs:
%   - pi_hat:    probability of choosing patch i from budget j
%   - N:        number of surveyed households in year N(t)

%% Global parameters
% This function is passed through a parfor loop.  Therefore, we cannot
% include global variables.  Instead, I put the variables that we need as
% an argument.  

%% Pre-assign variables
N = zeros(budget_l,1);        % number of people
N_patches = size(X,1);        % Number of patches
pi_hat = zeros(N_patches,1);  % probability of choosing patch i from budget j.
s{budget_l} = [];             % Polynomial basis function for pihat  
Q{budget_l} = [];             % Inner product of basis function: s.'*s 

%% Compute pi_hat.  
% We loop through each budget year
for jj = 1:budget_l
    % Number of people in year jj
    N(jj) = size(shares{jj},1);
    
    % Reset temporary variables
    patch = zeros(N(jj),budget_l);  % Patch that person n is on
    d = zeros(N(jj),1);             % Indicator to math person n to patch i
    
    % Implied quantity of each good demanded by person n:
    quantity = shares{jj}./kron(ones(N(jj),1),budgets(jj,:));
    % That is,q = shares/(p/e)=shares*e/p=(pq/e)/(p/e)=shares/budgets
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
    [~, d] = ismember(patch,X,'rows'); 

    % Construct basis function Z = [constant, income, income^2, ... , income^p]
    % We have already constructed the basis function for median income.
    s{jj} = ones(N(jj),1);      
    for pp = 1:poly_degree
        s{jj}  = [s{jj}, income{jj}.^pp];
    end
    
    % Q is the sum of outer products of the basis function s.  This is used
    % to weight pi_hat.
    Q{jj} = (s{jj}.'*s{jj});
    
    % Find all unique patches of X that people purchase from.  
    % We ignore patches not purchased from, since pi_hat is zero for those.
    uniq_patch = unique(d); 
    for ii = 1:size(uniq_patch)
        % Calculate pi_hat using formula in Yuichi's notes.  
        ind = d == uniq_patch(ii);
        pi_hat(uniq_patch(ii)) = med_income{jj}*(Q{jj}\sum(s{jj}(ind,:),1).');
    end
end

% Apply uniform G function to pi_hat, so it is between 0 and 1. 
pi_hat(pi_hat < 0) = 0;
pi_hat(pi_hat > 1) = 1;

% pi_hat does not sum to one within a year now.  Currently not normalized.
% See Joerg.

end
 






