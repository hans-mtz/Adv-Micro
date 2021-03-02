function  [KS_stat,KS_stat_b,KS_CV, KS_pval]  = RUM_21_KSTest(pi_hat,pi_b,X,A,RUMparams)
%% Code Description

%% Extract parameters
bootstrap_reps  = RUMparams.bootstrap_reps;   
seed            = RUMparams.seed;              
num_cores       = RUMparams.num_cores;  
ind_test        = RUMparams.ind_test;
N               = RUMparams.N;
I               = RUMparams.I;
J               = RUMparams.J;
ind_tau         = RUMparams.ind_tau;                                     

%% Preset output
KS_stat 	= zeros(1,1);
KS_stat_b   = zeros(1,bootstrap_reps);
KS_CV       = zeros(1,2);
KS_pval     = zeros(1,1);

%% Tuning parameter tau
% We set tau = sqrt(log(N)/N)) where N is the sample size.  This is in 
% contrast to the empirical example where we set N to be a function of the
% estimated variance and sample sizes across budgets.
% N_all is the total number of obesrvations over all budgets.
N_all   = N*J;
if ind_tau
    tau    = (log(N)/N)^0.5;
else
    tau = 0;
end

%% Jstat and bootstrap distribution
% Jstat is not tightened
[nuhat,Jstat] =  RUM_22_KSTeststat(A,pi_hat,N_all,0);

% Get tau-tightened nuhat
[nuhat_tight,~] =  RUM_22_KSTeststat(A,pi_hat,N_all,tau);

% Resample shares income to get non-parametric pihat.  
for bb = 1:bootstrap_reps
    % Set seed: Since we are doing this in parallel we need to set seed on 
    % each loop.  To ensure independence we set subseed to be B+bb rather 
    % than the seed (NB: B+bb is to avoid using same sub-seed as used
    % in computing tau).
    stream=RandStream('mlfg6331_64','Seed',seed);
    RandStream.setGlobalStream(stream);
    stream.Substream = bootstrap_reps + bb;

    % Recenter the tau-tightned nuhat
    pi_bs2 = pi_b(:,bb) - pi_hat + nuhat_tight;

    % Get the Jstatistic associated with the recenter the tau-tightned nuhat
    % Do not update tau.
    [~,Jstat_bs(bb,1)] =  RUM_22_KSTeststat(A,pi_bs2,N_all,tau);   
end

%% Critical value and Pr(Jstat == 0)
Jstat_bs = sortrows(Jstat_bs,1);
cv_95 = Jstat_bs( ceil(bootstrap_reps*0.95));
cv_99 = Jstat_bs( ceil(bootstrap_reps*0.99));
CV = [cv_95  cv_99];
if Jstat == 0
    prob = 1;
elseif isempty(min(Jstat_bs( Jstat_bs>=Jstat)))
    prob = 0;
else
    prob  = 1- (find(min(Jstat_bs( Jstat_bs>=Jstat)) == Jstat_bs,1)-1)/bootstrap_reps;
end


%% Output
KS_stat 	= Jstat;
KS_stat_b   = Jstat_bs.';
KS_CV       = CV;
KS_pval     = prob;

end