function [pi_hat,N] =RUM_73_EndogenousPihat(budgets,shares,income,med_income,Z,X,budget_l,poly_degree)
%% Code Description: Create pi_hat
% This module estimates the probability that patch patch i is chosen from
% budget set j, which we call pi_hat. 
% Inputs: 
%   - budgets:      median-income normalized budget lines
%   - shares:       share of expenditure on each composite good.
%   - income:       consumer income proxied with expenditure.
%   - med_income:   median income
%   - X:            Vector determining valid patches on each budget B_j
%   - budget_l:     budget length
%   - poly_degree:  polynomial degree
% Outputs:
%   - pi_hat:    probability of choosing patch i from budget j
%   - N:         number of surveyed households in year N(t)


%% Global parameters
% This function is passed through a parfor loop.  Therefore, we cannot
% include global variables.  Instead, I put the variables that we need as
% an argument.  

%% Pre-assign variables
N = zeros(budget_l,1);        % number of people
N_patches = size(X,1);        % Number of patches
pi_hat = zeros(N_patches,1);  % probability of choosing patch i from budget j.
r{budget_l} = [];             % Polynomial basis function for instrument
s{budget_l} = [];             % Polynomial basis function for pihat  
Q{budget_l} = [];             % s.'*s
epsilon_hat{budget_l} = [];   

%% Compute epsilon_hat
for jj = 1:budget_l
    % Sample size for period jj
    N(jj) = size(shares{jj},1);

    % Construct basis function for instrument 
    % r = [constant, Z, Z^2, ... , Z^p]
    r{jj} = ones(N(jj),1);      
    for pp = 1:poly_degree
        r{jj}  = [r{jj}, Z{jj}.^pp];
    end
    
    % Variance of basis function
    R_hat = (r{jj}.'*r{jj});

    % epsilon_tilde_{n(j)} is F(w|z) evalauted at w_{n(j)} and z_{n(j)},
    % where F(w|z) is given in endogeneity section.
    epsilon_tilde = zeros(N(jj),1);
    for nn = 1:N(jj)
        ind = income{jj}(nn,1) >= income{jj};
        epsilon_tilde(nn,1) =  r{jj}(nn,:)*(R_hat\(r{jj}.'*ind));
    end
    
    % Truncate epsilon_tilde to get epsilon_hat
    epsilon_hat{jj} = epsilon_tilde;
    v_N = 1/N(jj)^(1/3);
    iota_U =  (1 - epsilon_tilde + v_N).^2/(4*v_N);
    iota_L =  (epsilon_tilde + v_N).^2/(4*v_N); 
    
    % Truncation 1
    ind = epsilon_tilde > 1 + v_N;
    epsilon_hat{jj}(ind) = 1;
    
    % Truncation 2
    ind = epsilon_tilde > 1-v_N & epsilon_tilde <= 1 + v_N; 
    epsilon_hat{jj}(ind) =  1 - iota_U(ind);
    
    % Truncation 3
    ind = epsilon_tilde > -v_N & epsilon_tilde <= v_N; 
    epsilon_hat{jj}(ind) =  iota_L(ind);
    
    % Lower truncation
    ind = epsilon_tilde <  -v_N;
    epsilon_hat{jj}(ind) = 0;
end



%% Compute alpha.  
% We loop through each budget year
alpha{budget_l} = [];

for jj = 1:budget_l
    % Reset temporary variables
    patch = zeros(N(jj),budget_l);  % Patch that person n is on
    d = zeros(N(jj),1);             % Matches patch to X  
    alpha{jj} =  zeros(size(X,1),sum(1:poly_degree+1));
    
    % Implied quantity of each good demanded by person n:
    quantity = shares{jj}./kron(ones(N(jj),1),budgets(jj,:));
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
    [~, d] = ismember(patch,X,'rows'); 

    % Construct basis function 
    % s = constant, inc, eps, inc^2, inc*eps, eps^2, ...]
    % s = [constant, income, income^2, ... , income^p,epsilon_hat, epsilon_hat^2, ... , epsilon_hat^pp]
    s{jj} = zeros(N(jj),1);      
    for pp = 0:poly_degree
        for qq = 0:poly_degree
            if pp+qq <= poly_degree
                s{jj}  = [s{jj}, income{jj}.^pp.*epsilon_hat{jj}.^qq];
            end
        end
    end
    s{jj}(:,1) = []; 
    
    % S is the sum of outer products of the basis function s.  This is used
    % to weight pi_hat.
    S{jj} = (s{jj}.'*s{jj});
    
    % Find all unique patches of X that people purchase from.  
    % We ignore patches not purchased from, since pi_hat is zero for those.
    uniq_patch = unique(d); 
    for ii = 1:size(uniq_patch)
        % Calculate pi_hat using formula in Yuichi's notes.  
        ind = d == uniq_patch(ii);
        alpha{jj}(uniq_patch(ii),:) = (S{jj}\sum(s{jj}(ind,:),1).').';
    end
end

%% Compute pi_hat.  
pi_hat =  zeros(size(X,1),1);
D{budget_l} = [];
for jj = 1:budget_l
    % Construct basis function at median income and E(e):
    % Note that D  =(1, med(inc), int_0^1 e, med(inc)^2, med(inc)*int_0^1e
    %int_0^1 e^2, ...
    % Nb that int_0^1 e^p = 
    D{jj}  = zeros(1,1);  
    for pp = 0:poly_degree
        for qq = 0:poly_degree
            if pp + qq <= poly_degree
                D{jj} = [D{jj}, med_income{jj}.^pp/(qq+1)];
            end
        end
    end
    D{jj}(:,1) = [];
    pi_hat = pi_hat + (D{jj}*alpha{jj}.').';
end

% Apply uniform G function to pi_hat, so it is between 0 and 1. 
pi_hat(pi_hat < 0) = 0;
pi_hat(pi_hat > 1) = 1;

% pi_hat does not sum to one within a year now.  Currently not normalized.
% See Joerg.


end
 






