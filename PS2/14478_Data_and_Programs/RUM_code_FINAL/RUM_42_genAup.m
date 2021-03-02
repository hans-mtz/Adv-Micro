function A=RUM_42_genAup(column,level,counter,RP,A)
%% Code Description: Tree-crawling algorithm (2nd ... N steps)
% This is a sub-module for RUM_41_genA. 
% This function starts off with a vector defining a patch
% choice from the first budget.  Next, it adds a candidate patch vector
% from the second budget.  It runs the Floyd-Warshall algorithm to check if
% transitive closure holds.  If it holds, then the algorithm continues by
% calling itself and adding a new candidate patch.  It keeps going until
% either (1) transitive closure fails, in which case we ignore this
% candidate vector (and consequently all candidate vectors below this one
% on the tree) and try another candidate vector; or (2) we exhaust all the 
% budget years and have a valid pseudo-agent.
% Input:
%   - column:  candidate patches
%   - level:   how many patches we have added to column (i.e., keeps track
%              what budget-year we are in.
%   - counter: total number of patches for this level
%   - RP:      Revealed-preferred matrix used to check for cycles
%   - A:       matrix of pseudo-agents

if level ~= length(counter)
    for tt=1:counter(level) 
        posit = sum(counter(1:level-1)) + tt;
        columnnew = column;
        columnnew(posit)=true;
        ind_check = find(columnnew==1);
        AreCycles=RUM_43_genAfw( RP(ind_check,ind_check));
        if AreCycles==false %if no cycle detected, recursive call
            A=RUM_42_genAup(columnnew,level+1,counter,RP,A);
        end
    end
else % We are at the final level
    for tt=1:counter(level)
        posit = sum(counter(1:level-1)) + tt;
        columnnew = column;
        columnnew(posit)=true;
        ind_check = find(columnnew==1);
        AreCycles=RUM_43_genAfw( RP(ind_check,ind_check));
        if AreCycles==false %if no cycle detected, extend A
            A=[A,columnnew];
        end
    end
end
    
end   