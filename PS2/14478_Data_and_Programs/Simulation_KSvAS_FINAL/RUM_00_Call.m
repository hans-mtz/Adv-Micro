% Call the main file that executes the tests
% This file is a wrapper that sets the parameters of the model and calls
% RUM_00_Main_Simulations, which is the file that carries out the AS and KS tests.

% Clear data
clc
clear

% Turn warning off
warning('off','MATLAB:nargchk:deprecated')

% Inputs
bs_reps     = 999;        % Bootstrap reps
M           = 1000;       % Simulations
N           = 1000;       % Sample size
lambda      = [0:0.1:1];  % Weights on pi_inside and pi_outside
ind_test    = 3;          % Equal to 1 to run KS, 2 to run AS, 3 to run both
ind_tau     = 1;          % Set equal to 0 to force tau = 0 and kappa = inf
ind_33      = 0;          % Run DGP1-3 if set equal to 1, DGP4 else.

% Setup CVX
RUM_01_cvx

if ind_33
    % Load data
    load('Input\33_testdata_2.mat') %testdata_j is dgp for j'th panel of table

    % Update data
    A       = A33;
    B       = B33;
    Beq     = [];
    X       = X33;
    pi_lo   = pi_inside;
    pi_hi   = pi_outside_2;

    % Call tests for each lambda
    for ii =1:size(lambda,2)
        pi_true =  lambda(ii)*pi_lo + (1- lambda(ii))*pi_hi;

        RUM_00_Main_Simulation(A,X,B,Beq,pi_true,bs_reps,M,N,ind_test,lambda(ii),ind_tau)
    end
else
    % Load data
    load('Input\testdata_new.mat')

    % Update data
    A       = A;
    B       = B;
    Beq     = [];
    X       = X;

    % Call tests for each lambda
    RUM_00_Main_Simulation(A,X,B,Beq,pi_true,bs_reps,M,N,ind_test,nan,ind_tau)
end