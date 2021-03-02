function AreCycles=RUM_43_genAfw(F)
%% Code Description:  Floyd-Warshall Algorithm
% Takes an F matrix of integers (where n>0 in position i,j indicates 
% that jBi) and calculates the transitive closure using the 
% Floyd-Warshall algorithm.

n =length(F);
for kk=1:n
    F(bsxfun(@times,F(:,kk),F(kk,:)) == true) = true;
end

AreCycles=false;
if max(diag(F))==true % max(diag(F))==0 if there are no cycles
    AreCycles=true; 
end
end