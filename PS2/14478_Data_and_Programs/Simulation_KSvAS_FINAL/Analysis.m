% Post analysis file.
% After running the main simulation, call this file with the appropriate
% parameters defined below.  Summary statistics such as the average p-value
% of the bootstrap repetitions is reported.

clear
clc

cd Output\DGP1\N1000\

N_all       = 1000;
M_all       = 500;
B_all       = 499;
ind_tau     = 1;
lambda_all  = [0:0.1:1];
num_lambda  = size(lambda_all,2);


ASpval_mean = zeros(num_lambda,1);
ASstat_mean = zeros(num_lambda,1);
AScv_mean   = zeros(num_lambda,2);
ASpower_all = zeros(num_lambda,2);
keep_num    = zeros(num_lambda,1);
KSpval_mean = zeros(num_lambda,1);
KSstat_mean = zeros(num_lambda,1);
KScv_mean   = zeros(num_lambda,2);
KSpower_all = zeros(num_lambda,2);



for iter =1:num_lambda
    name = strcat('AS_and_KS_Test_B=',num2str(B_all),'_M=',num2str(M_all),'_N=',num2str(N_all),'_lambda=',num2str(lambda_all(iter)),'_indtau=',num2str(ind_tau),'.mat');
    load(name);
    ASpval_mean(iter,1)  = mean(AS_pval);
    ASstat_mean(iter,1)  = mean(AS_stat);
    AScv_mean(iter,:)    = mean(AS_CV);
    ASpower_all(iter,:)  = AS_power;
    keep_num(iter,:)     = mean(b_keep);
      
    KSpval_mean(iter,1)  = mean(KS_pval);
    KSstat_mean(iter,1)  = mean(KS_stat);
    KScv_mean(iter,:)    = mean(KS_CV);
    KSpower_all(iter,:)  = KS_power;
    
    clearvars -except ASpval_mean ASstat_mean AScv_mean ASpower_all keep_num KSpval_mean KSstat_mean KScv_mean KSpower_all N_all M_all B_all lambda_all iter ind_tau
end

KS_Table = [KSpower_all KSstat_mean KSpval_mean KScv_mean];
AS_Table  =[ASpower_all ASstat_mean ASpval_mean AScv_mean keep_num];
    
name = strcat('AS_and_KS_Test_B=',num2str(B_all),'_M=',num2str(M_all),'_N=',num2str(N_all),'.mat');
save(name,'ASpval_mean','ASstat_mean','AScv_mean','KSpval_mean','KSstat_mean','KScv_mean','KSpower_all','ASpower_all','keep_num','KS_Table','AS_Table')

cd ..
cd ..
cd ..