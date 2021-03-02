function [] = RUM_00_Main_Simulation(A,X,B,Beq,pi_true,bs_reps,M,N,ind_test,lambda,ind_tau)
%% Code Description: Simulation File
% Code by:
%	Yuichi Kitamura
%	Joerg Stoye
%	Matthew Thirkettle 
%
% This script computes the KS and AS statistics for (A,X) and (B,X) implied
% by the 3-by-3 empirical example.  We generate data using pi_true, where
% pi_true is the weighted sum of pi_inside and pi_outside.  pi_inside
% strongly satisfies the null Bpi_inside <= 0 .  On the other hand,
% pi_outside does not satisfy the null -- Bpi_outside > 0 for some
% components.  We map out power curves for the two statistics by varying
% pi_true from pi_inside to Bpi_outside.
%
% Reference paper: 
%	Non-parametric analysis of random utility models (RUM): testing.  
%	By Yuichi Kitamura and Joerg Stoye
%	Institute for Fiscal Studies, Departments of economics,
%   UCL Working Paper WP# CWP36/13 
%
% Last modified: October 24, 2017.

%% Start timer
tic  

%% Structure of parameters
RUMparams.ind_test          = []; 
RUMparams.bootstrap_reps    = [];
RUMparams.seed              = [];
RUMparams.num_cores         = [];
RUMparams.N                 = [];
RUMparams.I                 = [];
RUMparams.J                 = [];
RUMparams.ind_tau           = [];

%% Parameters
RUMparams.bootstrap_reps    = bs_reps;                                      % Number of bootstrap repetitions.   
RUMparams.N                 = N;                                            % Sample size                       
RUMparams.seed              = 10000*(RUMparams.bootstrap_reps/1000);        % Seed
RUMparams.ind_test          = ind_test;                                     % Indicator for test
[I,J]                       = size(X);                                      % Num patches and budgets
RUMparams.I                 = I;                                            % Total number of patches
RUMparams.J                 = J;                                            % Total number of budgets
RUMparams.ind_tau           = ind_tau;                                      % Indicator for tau = 0

%% Extract parameters           
bootstrap_reps              = RUMparams.bootstrap_reps;   
seed                        = RUMparams.seed;              
ind_test                    = RUMparams.ind_test;
N                           = RUMparams.N;
I                           = RUMparams.I;
J                           = RUMparams.J;
ind_tau                     = RUMparams.ind_tau;

%% Set seed
stream=RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream);
stream.Substream = 1;

%% Generate Data
% Determine number of patches in each budget
patch_size = zeros(J,1);
patch_ind  = zeros(J,1);
for jj = 1:J
    patch_size(jj,1) =  size(find(X(:,jj) == 0),1);
    patch_ind(find(X(:,jj) == 0),1) = jj;
end
patch_size = [0;patch_size];

% Preset data Choice.  Choice_{i,j,m} = d iff patch d is chosen from budget 
% j in Monte Carlo m.
Choice = zeros(N,J,M);
for jj = 1:J
    % Get pi_j
    c = sum(patch_size(1:jj));
    pi_j = pi_true(c+1:c+patch_size(jj+1));
    
    % Determine the CDF of pi_j
    F_j = [0;cumsum(pi_j)];
    
    % Uniform draw
    U = rand(N,M);
    
    % Map U to choices
    Choice_temp = zeros(N,M);
    for ii = 1:size(pi_j,1)
        Choice_temp( U > F_j(ii) & F_j(ii+1)>= U) = ii;
    end
    Choice(:,jj,:) = Choice_temp;
end

%% Estimate pi_hat and pi_hat_bs
% We use a simple frequency estimator.
pi_hat = zeros(I,M);
pi_b = zeros(I,M,bootstrap_reps);
for mm = 1:M
    % Frequency estimator for pi
    [pi_hat(:,mm)] = RUM_11_Pihat(Choice(:,:,mm), X);
    
    % Bootstrapped frequency estimator
    for bb = 1:bootstrap_reps
        stream=RandStream('mlfg6331_64','Seed',seed);
        RandStream.setGlobalStream(stream);
        stream.Substream = 1+bb;
        
        % Draw indicies
        ind = ceil(N*rand(N,1));
        
        % Set data
        Choice_b = Choice(ind,:,mm);
        
        % Compute pi_b
        [pi_b(:,mm,bb)] = RUM_11_Pihat(Choice_b, X);
    end
end


%% KS Test
% For each Monte Carlo, compute the KS test 
% Preset variables
if ind_test == 1 || ind_test == 3
    KS_stat     = zeros(M,1);
    KS_stat_b   = zeros(M,bootstrap_reps);
    KS_CV       = zeros(M,2);
    KS_pval     = zeros(M,1);
    
    % Compute KS test statistic
    parfor mm = 1:M
        [KS_stat(mm,1),KS_stat_b(mm,:),KS_CV(mm,:), KS_pval(mm,1)]  = RUM_21_KSTest(pi_hat(:,mm),pi_b(:,mm,:),X,A,RUMparams);
    end
end

%% AS Test
if ind_test == 2 || ind_test == 3
    % For each Monte Carlo, compute the AS test 
    % Preset variables
    AS_stat     = zeros(M,1);
    AS_stat_b   = zeros(M,bootstrap_reps);
    AS_CV       = zeros(M,2);
    AS_pval     = zeros(M,1);
    b_drop      = zeros(M,1);
    % Compute AS test statistic
    parfor mm = 1:M
       [AS_stat(mm,1),AS_stat_b(mm,:),AS_CV(mm,:), AS_pval(mm,1), b_keep(mm,1)]  = RUM_31_ASTest(pi_hat(:,mm),pi_b(:,mm,:),X,B,Beq,RUMparams);
    end
end


%% Power
% Determine power of the tests KS and AS
if ind_test == 1 || ind_test == 3
    % Power is equal to average number of times we reject H0 at alpha =
    % 0.05 or 0.01.  That is, it is the average number of times KS_stat
    % exceeds KS_CV.
    KS_power = zeros(1,2);
    KS_power(1,1) = size(find(KS_stat > KS_CV(:,1)),1)/M;
    KS_power(1,2) = size(find(KS_stat > KS_CV(:,2)),1)/M;
end
if ind_test == 2 || ind_test == 3
    AS_power = zeros(1,2);
    AS_power(1,1) = size(find(AS_stat > AS_CV(:,1)),1)/M;
    AS_power(1,2) = size(find(AS_stat > AS_CV(:,2)),1)/M;
end

%% Save output
Time_taken = toc;
c = clock;
datetime = strcat(num2str(c(1)),num2str(c(2)),num2str(c(3)),'_',num2str(c(4)),num2str(c(5)));
name = strcat('Output/AS_and_KS_Test_B=',num2str(bootstrap_reps),'_M=',num2str(M),'_N=',num2str(N),'_lambda=',num2str(lambda),'_indtau=',num2str(ind_tau),'.mat');
save(name); 
end




