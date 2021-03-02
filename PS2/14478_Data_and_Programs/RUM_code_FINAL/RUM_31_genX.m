function X=RUM_31_genX(budgets)
%% Code Description: Generate X
% This code generates a list of existing "patches" X from a set of budgets.
% We build a matrix X of all possible patches.  We delete patches that are
% not valid, i.e., not in the positive quadrant.
% Input: budgets = [budget_l,classes] matrix of normalized budgets.
% Output: X = a [?,budget_l] matrix listing the existing budgets as 
% "over-under-on" strings.

% Globals
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores

%% Generate all potential possible combinations of X
% Pre-allocate X to be zero matrix (int8 as elements can be in {-1,0,1}
[J,dimx]=size(budgets);
X=int8(zeros(J*2^(J-1),J));

% Generate a [?,J-1] matrix of (-1,1) combinations
referencematrix=int8(zeros(2^(J-1),J-1));
for tt=1:J-1
    block=int8([-ones(2^(tt-1),1);ones(2^(tt-1),1)]);
    for count=1:2^tt:2^(J-1)-2^tt+1
        referencematrix(count:count+2^tt-1,tt)=block;
    end
end

% Generate all logically possible (0,-1,1)-strings
for tt=1:J
    block=int8(zeros(2^(J-1),J));
    if tt>1
        block(:,1:tt-1)=referencematrix(:,1:tt-1);
    end
    if tt<J
        block(:,tt+1:J)=referencematrix(:,tt:J-1);
    end
    X(2^(J-1)*(tt-1)+1:2^(J-1)*tt,:)=block;
end

%% Test and Trim X 
% Number of rows to check 
I = J*2^(J-1);

% Number of loops
L = ceil(I/num_cores);

Xdrop{num_cores} = [];
parfor cc = 1:num_cores
    Xdrop2 = -999*ones(L,1);
    for ll = 1:L
        % Row to check
        ii = (cc-1)*L + ll;
        if ii > I
            continue;
        end
        Xdrop2(ll,1) =  RUM_32_genXdroprow(budgets,X,dimx,ii,J);
    end
    Xdrop{cc} = Xdrop2.';
end
xdrop = [Xdrop{:}].';
xdrop(xdrop == -999,:) = [];

%Drop invalid patches (rows of x)
X(xdrop(:,1) == 1,:) = [];

    
end
        





        