function [AS_stat,AS_stat_b,AS_CV, AS_pval,num_rows_keep] =  RUM_31_ASTest(pi_hat,pi_b,X,B,Beq,RUMparams);
%% Code Description
% Assume that we know an M*I matrix B is stochastic rationaliability of pi
% is equivalent to Bpi <=0.  Also, we have an estimator pi_hat whose
% bootstrap analog pi_hat* is defined iid resampling of the data with
% replacement.  The distribution of pi_hat* is simulated through drawing B
% = 1001 bootstrap realizations pi_bs_b.
%
% In this function we construct a test statistic and bootstrap
% distribution for testing B*pi <= 0.
%
% Input:
%   - B and Beq     M-by-I and Meq-by-I matrices representing inequality
%                   and equality tests.  Under the null, B*x <=0 and
%                   Beq*x=0.
%   - pi_hat.       This is I*1.
%   - pi_bs.        Bootstrap distribution for pi_hat*
%   - RUMparams.    Set of parameter options
%   - shares.       Share of person i's income spent on good k
%
% Outputs:
%   - T_N.          Test statistic
%   - T_bs.         Bootstrap distribution for T_N.

%% Extract parameters           
bootstrap_reps  = RUMparams.bootstrap_reps;   
seed            = RUMparams.seed;              
num_cores       = RUMparams.num_cores;  
ind_test        = RUMparams.ind_test;
N               = RUMparams.N;
I               = RUMparams.I;
J               = RUMparams.J;
ind_tau         = RUMparams.ind_tau;

%% Step 0
% Sample size
N_all      = N*J;

% Size B
b_size = size(B,1);

% Determine the number of patches for each budget j
% Determined by X.
patch_size = zeros(J,1);
patch_ind  = zeros(J,1);
for jj = 1:J
    patch_size(jj,1) =  size(find(X(:,jj) == 0),1);
    patch_ind(find(X(:,jj) == 0),1) = jj;
end

%% Step 1: Omega_hat
% Compute estimated variances of moment conditions
% For each budget-year j = 1, ... , J we compute block j in the block diagonal
% matrix Sigma_hat, where
% Sigma_hat_j = diag(pi_hat_j) - pi_hat_j*pi_hat_j'.
Sigma_hat = zeros(I,I);
pi_j = pi_hat(1:patch_size(1));
Sigma_hat(1:patch_size(1),1:patch_size(1)) = diag(pi_j) - pi_j*pi_j.';
for jj = 2:J
    c = sum(patch_size(1:jj-1));
    pi_j = pi_hat(c+1:c+patch_size(jj));
    Sigma_hat(c+1:c+patch_size(jj),c+1:c+patch_size(jj)) ...
        = diag(pi_j) - pi_j*pi_j.';
end

% Estimated variance is the usual sandwhich formula Omega = B*Sigma_hat*B'
% Get diagonal elements
if ~isempty(B) && isempty(Beq)
    Omega       = B*Sigma_hat*B.';
    Omega_diag      = Omega(find(eye(size(Omega))));
    Omega_eq_diag   = [];
elseif isempty(B) && ~isempty(Beq)
    Omega_eq        =  Beq*Sigma_hat*Beq.';
    Omega_eq_diag   = Omega_eq(find(eye(size(Omega_eq))));
    Omega_diag      = [];
else
    Omega       = B*Sigma_hat*B.';
    Omega_diag      = Omega(find(eye(size(Omega))));
    Omega_eq        =  Beq*Sigma_hat*Beq.';
    Omega_eq_diag   = Omega_eq(find(eye(size(Omega_eq))));
end

%% Step 2: 
% Drop inequalites that are not close to binding.
% We drop rows that violate the condition
% sqrt(N)*B^m*pi_hat*Omega_hat_{mm}^{-1/2} >= -kappa_N.

% First, check that diagonal elements of Omega_diag are positive
if ~isempty(find(Omega_diag<=0, 1)) || ~isempty(find(Omega_eq_diag<=0, 1))
   error('Diagonal element of Omega is non-positive')
end

if ~isempty(B)
    % Determine if there exists any rows that violate the GMS condition
    if ind_tau
        kappa_N = sqrt(log(N));
    else
        kappa_N = 0;
    end

    val_tilde       = sqrt(N)*(B*pi_hat).*(Omega_diag.^(-0.5)) + kappa_N;
    ind_tilde       = find(val_tilde >= 0);
    if isempty(ind_tilde) && isempty(Beq)
        AS_stat 	= 0;
        AS_stat_b   = zeros(1,bootstrap_reps);
        AS_CV       = zeros(1,2);
        AS_pval     = 1;
        num_rows_keep = 0;
        warning('There does not exist a row of B_tilde satisfying the kappa_n condition and Beq is empty')
        return
    end
    B_tilde         = zeros(size(ind_tilde,1),I);
    B_tilde         = B(ind_tilde,:);

    % Update B and Omega
    B               = B_tilde;
    Omega           = B*Sigma_hat*B.';
    Omega_diag      = Omega(find(eye(size(Omega))));
end

%% Step 3: Test statistic
% The MMM test statistic is 
% T_N  = N*(sum_{m} max(0,B^m*pi_hat)^2*Omega_{mm}^(-1/2)
%      + N*(sum_{m} (Beq^m*pi_hat)^2*Omegaeq_{mm}^(-1/2)
if ~isempty(B) && isempty(Beq)
    B_pi    = max(B*pi_hat,0).^2;
    T_N     = N_all*sum((B_pi).*(Omega_diag.^(-0.5)));
elseif isempty(B) && ~isempty(Beq)
    Beq_pi  = (Beq*pi_hat).^2;
    T_N     = N_all*sum((Beq_pi).*(Omega_eq_diag.^(-0.5)));
else
    B_pi    = max(B*pi_hat,0).^2;
    Beq_pi  = (Beq*pi_hat).^2;
    T_N     = N_all*sum((B_pi).*(Omega_diag.^(-0.5)))  +  N_all*sum((Beq_pi).*(Omega_eq_diag.^(-0.5)));
end

T_N = round(T_N,6,'decimals');

%% Step 4: Test statistic
% Bootstrap distribution of T_N
T_bs = zeros(bootstrap_reps,1);
for bb=1:bootstrap_reps
    pi_bs       = pi_b(:,bb);
    pi_r       = pi_b(:,bb) - pi_hat;
    
   
    % Calculate Omega_b
    Sigma_b = zeros(I,I);
    pi_j = pi_bs(1:patch_size(1));
    
    Sigma_b(1:patch_size(1),1:patch_size(1)) = diag(pi_j) - pi_j*pi_j.';
    for jj = 2:J
        c = sum(patch_size(1:jj-1));
        pi_j = pi_bs(c+1:c+patch_size(jj));
        Sigma_b(c+1:c+patch_size(jj),c+1:c+patch_size(jj)) ...
        = diag(pi_j) - pi_j*pi_j.';
    end

    if ~isempty(B) && isempty(Beq)
        Omega_b         = B*Sigma_b*B.';
        Omega_diag_b    = Omega_b(find(eye(size(Omega_b))));
        Omega_eq_diag_b = [];
    elseif isempty(B) && ~isempty(Beq)
        Omega_eq_b      =  Beq*Sigma_b*Beq.';
        Omega_eq_diag_b = Omega_eq_b(find(eye(size(Omega_eq_b))));
        Omega_diag_b    = [];
    else
        Omega_b         = B*Sigma_b*B.';
        Omega_diag_b    = Omega_b(find(eye(size(Omega_b))));
        Omega_eq_b      = Beq*Sigma_b*Beq.';
        Omega_eq_diag_b = Omega_eq_b(find(eye(size(Omega_eq_b))));
    end
    
    if ~isempty(find(Omega_diag_b<=0, 1)) || ~isempty(find(Omega_eq_diag_b<=0, 1))
        error('Diagonal element of Omega_b is non-positive')
    end

    if ~isempty(B) && isempty(Beq)
        B_pi_b    = max(B*pi_r,0).^2;
        T_bs(bb,1)= N_all*sum((B_pi_b).*(Omega_diag_b.^(-0.5)));
    elseif isempty(B) && ~isempty(Beq)
        Beq_pi_b  = (Beq*pi_r).^2;   
        T_bs(bb,1)=  N_all*sum((Beq_pi_b).*(Omega_eq_diag_b.^(-0.5)));
    else
        B_pi_b    = max(B*pi_r,0).^2;
        Beq_pi_b  = (Beq*pi_r).^2;   
        T_bs(bb,1)= N_all*sum((B_pi_b).*(Omega_diag_b.^(-0.5)))  +  N_all*sum((Beq_pi_b).*(Omega_eq_diag_b.^(-0.5)));
    end
end

T_bs = round(T_bs,6,'decimals');

%% Results
T_bs  = sort(T_bs);
cv_95 = T_bs( ceil(bootstrap_reps*0.95));
cv_99 = T_bs( ceil(bootstrap_reps*0.99));
CV = [cv_95  cv_99];
if T_N == 0
    prob = 1;
elseif isempty(min(T_bs( T_bs>=T_N)))
    prob = 0;
else
    prob  = 1- (find(min(T_bs( T_bs>=T_N)) == T_bs,1)-1)/bootstrap_reps;
end

%% Output
AS_stat 	= T_N;
AS_stat_b   = T_bs.';
AS_CV       = CV;
AS_pval     = prob;
num_rows_keep = size(B,1);
end