function [B_out,Beq_out] = RUM_genB(A,Rpath)
%% Code Description
% This function computes H representation of A using an rcdd.  This is
% called via R.
%
% The output is a Mineq-by-I matrix B_out indicating inequalities
% B_out*x<=0 and a Meq-by-I matrix Beq_out indicating equalities
% Beq_out*x=0.  Throughout this code, in constrast, B is the concatination of the
% equality and inequality constraints.  Beq define which rows are equality
% and which rows are inequality.
%
% INPUT:
%   A       Matrix A defining psuedo-agents (vertices)
%   Rpath   Directory to R.  Format: '"C:\Program Files\R\R-3.4.2\bin\x64"'
%
% OUTPUT:
%   B_out       Inequality constraints defining cone
%   Beq_out     Equality constraints defining cone



%% Save tempoary file
save('A_Vrepresentation.mat')

%% Call R file
path1 = getenv('PATH');
if ~contains(path1,Rpath)
    path1 = strcat(path1,Rpath,';');
    setenv('PATH', path1)
end
system('R CMD BATCH Rum_genB_VtoH.R');


%% Load and delete results.
delete('A_Vrepresentation.mat')
load('B_Hrepresentation.mat')
delete('B_Hrepresentation.mat')
delete('.RData')

%% Reduce matrix B
% Step 1) 
% Remove rows in B that correspond to non-negativity constraints x_i >= 0.
[M,I]       = size(B);
ind_drop    = zeros(M,1);
for mm = 1:M
    if sum(B(mm,:) == 0) >= I-1 && Beq(mm,1) == 0
        ind_drop(mm,1) = 1;
    end
end
B(ind_drop == 1,:) = [];
Beq(ind_drop == 1,:) = []; 

% Set B = [C;D] where C correspond to equality constraints and D
% corresponds to inequality constraints.    

% Step 2)
% Put C in canonical form using rref.
C       = B(Beq == 1,:);
D       = B(Beq == 0,:);
C       = rref(C);

% Step 3) 
% Put D in canonical form using linear combinations of C.
Mc      = size(C,1);
D_ind   = zeros(size(D,1),1); 
for ii = 1:I
    % Precheck that C is a pivotal column
    if size(find(C(:,ii)==0),1) ~= Mc - 1
        break
    end
    
    % Sort rows of D 
    ind2    = find(D_ind ~= 0);
    ind1    = find(D(:,ii) ~=0 & D_ind == 0);
    ind0    = find(D(:,ii) ==0 & D_ind == 0);
    D       = D([ind2;ind1;ind0],:);
    D_ind   = D_ind([ind2;ind1;ind0],1);
    ind1    = find(D(:,ii) ~=0 & D_ind == 0);
    
    % If D is a vector of zeros, then we can skip this column
    if isempty(ind1)
        continue
    end
    
    % Find row in C == 1
    mm = find(C(:,ii)==1);
    
    % Check to make sure row mm is pivotal
    ind = find(C(mm,:) ~= 0 );
    if ind(1) ~= ii
        % Row mm is not pivotal and D is not a vector of zeros. 
        % Flag rows of D that cannot be updated
        D_ind(D(:,ii) ~=0,1) = 1;
        continue
    end
    
    % Reduce rows in D
    cc = C(mm,ii);
    for jj = 1:size(ind1,1)
        dd = D(ind1(jj),ii);
        D(ind1(jj),:)  = D(ind1(jj),:) - (dd/cc)*C(mm,:);
    end
end
B   =[ C;D];
Beq = [ones(size(C,1),1);zeros(size(D,1),1)];

% Step 4)
% Remove rows of B that correspond to x_i >= 0.
[M,~]       = size(B);
ind_drop    = zeros(M,1);
for mm = 1:M
    if sum(B(mm,:) == 0) >= I-1 && Beq(mm,1) == 0
        ind_drop(mm,1) = 1;
    end
end
B(ind_drop == 1,:) = [];
Beq(ind_drop == 1,:) = []; 

% Output B_out and Beq_out
B_out   = B(Beq == 0,:);
Beq_out = B(Beq == 1,:);

end