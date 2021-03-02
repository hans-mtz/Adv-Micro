function [pi_hat, Jstat, Jstat_bs, CV, pval] =  RUM_61_KernelStatistics(budgets,shares,income,X,A) 
%% Code Description: Test statistic for H_0: Households satisfy SARP
% This module computes pihat (probability of choosing patch i from budget j),
% the Jstat, and the boostrapped Jstat distrubiton
% using Kernel smoothing (Silverman's rule of thumb). 
% First, compute pihat using Kernel Smoothing.
% Second, get untightened Jstat and tightened nu_hat using tau.
% Third, re-draw the data and calculate a new bootstrapped pihat.  
% Recenter the tightened nuhat using (pihat_bs - pihat).  
% Last, get critical values and probabilities for the untightened Jstat
% Input:
%   - budgets: median-normalized budget constraints
%   - shares:  Expenditure share for each person n
%   - income:  Total expenditure on non-durables for each person n
%   - X:       Matrix defining patches
%   - A:       Matrix of GARP consistent pseudo-agents 
% Output:
%   - pi_hat:   Probability of purchasing from patch i on budget j
%   - Jstat:    The tau-tightened J statistic
%   - Jstat_bs: Bootstrapped distribution of Jstat
%   - CV:       CV(1) is the 95 percentile critical value and CV(2) is the
%               99 percentile
%  - prob:      Probability that the tau-tightened J statistic is zero

%% Global parameters
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores tau_val


%% Pre-assign variables
Jstat_bs  = zeros(bootstrap_reps,1);

%% Log income
for jj = 1:budget_l
   income{jj} = log(income{jj}); 
end


%% Median income
% Calculate median income.  
med_income = zeros(budget_l,1);  
for jj = 1:budget_l
    med_income(jj,1) = median(income{jj});
end

%% Pihat
% We calculate pi_hat, which is the probability of a household with
% median income purchasing patch i from budget j.
[pi_hat,N_j,h_N] =RUM_63_KernelPihat(budgets,shares,income,med_income,X,budget_l);

%% Jstat and bootstrap distribution
% Tuning parameter
Nh = min(N_j.*h_N);
if tau_val == 1
    tau = (log(Nh)/Nh)^0.5; % Tau is fixed for bootstrap
else 
    tau = 0;
end

% Jstat is not tightened
[~,Jstat] =  RUM_64_KernelTeststat(A,pi_hat,N_j,0);

% Get tau-tightened nuhat
[nuhat_tight,~] =  RUM_64_KernelTeststat(A,pi_hat,N_j,tau);

% Resample shares income to get non-parametric pihat.  
% Get new n_k and compute re-centered Jstat bootstrap
J = ceil(bootstrap_reps/num_cores);
I = num_cores;
B = bootstrap_reps;
Jstat_bs_temp{I} = [];
Jstat_bs  = zeros(J*I,1);
parfor ii = 1:I
    Jstat_bs_temp2 = zeros(J,1);
    for jj = 1:J
        % Effective bootstrap iteration
        bb = (ii-1)*J + jj

        % If we have reached total bootstrap_reps, continue
        if bb > B
            continue;
        end
        % Set seed: Since we are doing this in parallel we need to set seed on 
        % each loop.  To ensure independence we set subseed to be bb rather 
        % than the seed.
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = bb;

        % Resample shares and income
        % We need to run this sub-function because MatLab cannot handle
        % structures in parallel.
        [shares_bs,income_bs] = RUM_62_KernelRandomize(shares,income,N_j,budget_l);

        % Get new pi_hat based on resampled shares and income, but do not 
        % update median income.
        pi_bs2 =RUM_63_KernelPihat(budgets,shares_bs,income_bs,med_income,X,budget_l);

        % Recenter the tau-tightned nuhat
        pi_bs2 = pi_bs2 - pi_hat + nuhat_tight;

        % Get the Jstatistic associated with the recenter the tau-tightned nuhat
        % Do not update tau.
        [~,Jstat_bs_temp2(jj,1)] =  RUM_64_KernelTeststat(A,pi_bs2,N_j,tau);
    end
    Jstat_bs_temp{ii} = Jstat_bs_temp2;
end
for ii=1:I
    Jstat_bs( (ii-1)*J+1:ii*J,1) = Jstat_bs_temp{ii};
end
if B ~= I*J
    Jstat_bs = Jstat_bs(1:B,1);
end

%% Critical value and Pr(Jstat == 0)
Jstat_bs = sortrows(Jstat_bs,1);
cv_95 = Jstat_bs( ceil(bootstrap_reps*0.95));
cv_99 = Jstat_bs( ceil(bootstrap_reps*0.99));
CV = [cv_95  cv_99];
if Jstat == 0
    pval = 1;
elseif isempty(min(Jstat_bs( Jstat_bs>=Jstat)))
    pval = 0;
else
    pval  = 1- (find(min(Jstat_bs( Jstat_bs>=Jstat)) == Jstat_bs,1)-1)/bootstrap_reps;
end


end