
for iter = 1:budget_n
    J = Jstat(iter,1);
    J_bs = Jstat_bs{iter};
    
    J = round(J,6,'decimal');
    J_bs = round(J_bs,6,'decimal');
    
    if J == 0
        pval(iter,1) = 1;
    elseif isempty(min(J_bs( J_bs>=J)))
        pval(iter,1) = 0;
    else
        pval(iter,1)  = 1- (find(min(J_bs( J_bs>=J)) == J_bs,1)-1)/bootstrap_reps;
    end
    
    Jstat(iter,1) = J;
    
end