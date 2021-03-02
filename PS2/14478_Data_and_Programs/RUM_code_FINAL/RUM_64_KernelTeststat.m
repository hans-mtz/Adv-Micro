function [nuhat,fval] =  RUM_64_KernelTeststat(A,pi_hat,N_j,tau,Omega) 
%% Code Description
% nuhat solves min (A*nu - pihat)'Omega(A*nu - pihat) subject to nu >= 0. 
% Input:
%   - A is the agent-type matrix that determines which patch a particular
%       type of agent picks from budget j.
%   - pihat is computed using data on patch choices
%   - N is a Jx1 vector with the number of people sampled from budget
%       j=1,...,J
%   - tau is the tuning parameter.  If not specified, it is set to 0.
%   -Omega is a consistent estimator for the asymptotic variance.  Default
%   is to set it equal to I.

%% Global parameters
% This function is passed through a parfor loop.  Therefore, we cannot
% include global variables.  Instead, I put the variables that we need as
% an argument.  

warning('off','MATLAB:nargchk:deprecated')
%% Setup
[I,H]=size(A);
if nargin<5
    Omega=eye(I);
end
if nargin<4
    tau=0;
end

% Check to see if Omega is PD and symmetric
[~,p] = chol(Omega);
if p > 0
    error('Omega is not positive definite')
end
if Omega ~= Omega.'
    error('Omega is not symmetric')
end
% N = sum N_j
N = sum(N_j);

% Lower-bound for v
low=tau.*ones(H,1)/H;

%% Find nuhat
cvx_begin quiet
    variable v(H)
    minimize( norm(chol(Omega)*(A*v-pi_hat),2)) % You have to minimize a norm in CVX, not a quadratic form.
    subject to
       v >= low;
cvx_end

%% Output
nuhat =A*v;
fval  =N*(A*v-pi_hat).'*Omega*(A*v-pi_hat);
fval  = round(fval,6,'decimals');

end
