function [shares_bs,income_bs] = RUM_62_KernelRandomize(shares,income,N,budget_l)
% This function generates a new shares and income dataset for each year
% using U[0,1] sampling with replacement.
% We need this subfunction due to complications with parfor - MatLab cannot
% handle structures in the parfor loop.

for ii = 1:budget_l
    shares_bs{ii} = zeros(N(ii),1);
    income_bs{ii} = zeros(N(ii),1);
    ind = ceil(N(ii)*rand(N(ii),1));
    shares_bs{ii} = shares{ii}(ind,:);
    income_bs{ii} = income{ii}(ind,:);
end

end