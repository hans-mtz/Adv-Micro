## All functions needed for computations

#This function computes the smallest number of observations per menu (needed to τ_n)
## X=is a dataset with columns menu,choice, and rows individuals
## There are 32 menus in our setup.
function nintaun(X)
    dM=maximum(X[:,1])
    Ntau=10000.0
    for i=2:dM # First menu is empty
      Ntau=minimum([Ntau,sum(X[:,1].==i)])
    end
    return Ntau
end

# Finding non-redundant indices
function gindexes(Menus,dYu)
    dM=length(Menus)
    MM=zeros(dM,dYu)
    if dYu==6
        #For RUM we add the default to every menu
        Menus=[vcat(1,Menus[i].+1) for i in eachindex(Menus)]
        Menus[1]=[]
    end
    for i in 1:dM
        MM[i,Menus[i][1:end-1]].=1.0 # The last option in the menu is dropped
    end
    # Nonzero and linearly independent frequencies in calibrated probabilities
    gindex=findall(MM.==1.0)
    return gindex
end

# Given data computes empirical frequence for all menus and options
function frequences(data)
    dM=maximum(data[:,1]); dY=maximum(data[:,2])-minimum(data[:,2])+1;
    F=zeros(Float64,dM,dY)
    for i in 1:size(data,1)
        F[data[i,1],data[i,2]+1]= F[data[i,1],data[i,2]+1]+1.0
    end
    # Computing sample frequences
    P=zeros(Float64,dM,dY)
    P[1,1]=1.0
    P[2:end,:]=F[2:end,:]./sum(F[2:end,:],dims=2)
    # Throwing away P(default) since it is equal to 1-sum(P[1:5])
    P=P[:,2:end]
    return P
end

# This function translates integer menu identifiers to actual menus
function int2menu(number, menus=Menus)
    return menus[number]
end

# This function translates actual menus to integer menu identifiers
function menu2int(vector::Array{Int64,1})
    p=true; i=1;
    while p
        if int2menu(i)==vector
            p=false
        end
        i=i+1
    end
    return i-1
end

# This function computes all subsets of a given set
# Uses Combinatorics
function subsetsint(intset,menus=Menus)
    Temp=collect(powerset(int2menu(intset, menus)))
    return [menu2int(Temp[i]) for i in 1:length(Temp)]
end

# Compute m_A(D) given the matrix of frequences ps
function m(D,A,ps)
    if A==1
        return "error"
    end
    if ~issubset(D,subsetsint(A))
        return 0.0
    else
        etavecT=etavec_con(ps) #Computing the attention index
        if model=="LA" # Logit atention model
            return etavecT[D]/sum(etavecT[subsetsint(A)])
        elseif model=="RCG"
            beta=0.0   # Random categorization model
            for i in 1:length(etavecT)
                if intersect(int2menu(i),int2menu(A))==int2menu(D)
                    beta=beta+etavecT[i]
                end
            end
            return beta
        end
    end
end

# Compute the atention index eta(D) for all D given the matrix of frequences ps
# This function uses the direct formula from the paper
function etavec(ps)
    # Compting P(default)
    po=1.0 .- sum(ps,dims=2)
    DD=subsetsint(32) # All subsets
    etavec=zeros(Float64,length(DD))
    for i=1:length(DD)
        BB=subsetsint(DD[i]) # Subsets of a given set
        Dt=int2menu(DD[i])
        betatemp=0.0
        for j=1:length(BB)
            betatemp=betatemp+elmobinv(Dt,BB[j],po) #adding up elements of inversion
        end
        etavec[i]=betatemp
    end
    return etavec
end

# This function computes the elements of the summation in the defintion if
# the atention index eta is directly computed
function elmobinv(Dt,BB,po)
    if model=="LA"
        return (-1.0)^(length(setdiff(Dt,int2menu(BB))))*po[32]/po[BB]
    elseif model=="RCG"
        return (-1.0)^(length(setdiff(Dt,int2menu(BB))))*po[menu2int(setdiff(vec(1:5),int2menu(BB)))]
    end
end


# Compute the constrained atention index eta(D) for all D given the matrix of frequences ps
# This function restricts etas to be probabilities
function etavec_con(ps)
    # Computing P(default)
    po=1.0 .- sum(ps,dims=2)
    DD=subsetsint(32) # All subsets
    etamin=Model(Gurobi.Optimizer)
    #set_optimizer_attribute(etamin,"outlev",0)
    @variable(etamin, etaparam[1:length(DD)]>=0)
    @constraint(etamin, addm, sum(etaparam[t] for t in 1:length(DD))==1)
    if model=="LA"
        ## p(o,X)/p(o,A)
        @objective(etamin,Min,sum((po[32]/po[DD[i]]-sum(etaparam[l] for l in subsetsint(DD[i])))^2 for i in 1:length(DD)))
    elseif model=="RCG"
        ## p(o,X-A)
        @objective(etamin,Min,sum((po[menu2int(setdiff(vec(1:5),int2menu(DD[i])))]-sum(etaparam[l] for l in subsetsint(DD[i])))^2 for i in 1:length(DD)))
    end

    JuMP.optimize!(etamin)

    return value.(etaparam)
end

# Computing sum_C m_A(C)p_{\pi}(a,C)
# To compute m the actual data p is used. p_\pi is captured by pp
function expectedm(a,A,p,pp)
    DD=subsetsint(A)
    em=0.0
    for i in 1:length(DD)-1 # -1 to takes into account that summation is taken over strict subsets.
        em=em+m(DD[i],A,p)*pp[DD[i],a]
    end
    return em
end

# Computing the calibrated p_\pi from the matrix of frequencies using the direct formula
function pipie(x)
    dM,dYu=size(x)
    pp=zeros(dM,dYu)
    for i in 2:dM # i=1 is an empty menu
        if length(int2menu(i))==1 #Finding singleton menus
            pp[i,int2menu(i)[1]]=1.0
        else
        for k in 1:dYu
            pp[i,k]=(x[i,k]-expectedm(k,i,x,pp))./m(i,i,x) #recursive formula
        end
        end
    end
    ## normalization to add up to 1
    pp[2:end,:]=pp[2:end,:]./sum(pp[2:end,:],dims=2)
    return pp
end

# Computing the constrained calibrated p_\pi from the matrix of frequencies
# This function is constrained to return probabilities
function pipie_cons(x)
    dM,dYu=size(x)
    MM=[ m(D,A,x) for D in 1:dM, A in 2:dM] # Matrix of consideration probabilites
    ConsB=(x.>0.0) # Matrix of 0/1. ConsB=0 if P(a in A)=0 and =1 otherwise
    pipiemin=Model(Gurobi.Optimizer)
    #set_optimizer_attribute(pipiemin,"outlev",0)
    @variable(pipiemin, pipieparam[1:dM,1:dYu]>=0)
    @constraint(pipiemin, sump[l=2:dM], sum(pipieparam[l,t] for t in 1:dYu)==1)
    @constraint(pipiemin, [l=1:dM,t=1:dYu], pipieparam[l,t]<=ConsB[l,t])

    @objective(pipiemin,Min, sum([(x[A,a]-sum(MM[D,A] * pipieparam[D,a] for D in 1:dM))^2 for A in 1:dM-1, a in 1:dYu]))
    JuMP.optimize!(pipiemin)

    return value.(pipieparam)
end


# Computing G
# U is the set of preferences
# M is the set of menus
# gindexsh are coordinates of nonzero linearly independent p_pi
function matrixcons(gindexsh, M, U)
    dYu=length(U[1])
    dYm=length(M)
    if dYu==6
        M=[vcat(1,M[i].+1) for i in eachindex(M)]
        M[1]=[]
    end
    d=length(U)
    d2=length(gindexsh)
    B=zeros(d2,d)
    m1=1
    for j in 1:d
        pp=zeros(dYm,dYu)
        for i in eachindex(M)
            if length(M[i])>0 # Skipping empty menu
                for k in 1:dYu
                    if k==M[i][argmax(U[j][M[i]])]
                        pp[i,k]=1.0
                    end
                end
            end
        end
        B[:,m1]=pp[gindexsh] #picking only relevant elements, indexes that are not always zero
        m1=m1+1
    end
    if dYu==6 # If model is RUM then return B
        return B
    else #otherwise return
        return [[B zeros(size(B,1),dYm)];[zeros(dYm,size(B,2),) I]]
    end
end

# This function computes the vector if linearly independent nozero elements of
# p_\pi and m (if needed)
function estimateg(X,gindex)
  Y=frequences(X)
  if length(gindex)==80 # If model is RUM
       return [1.0 .- sum(Y,dims=2) Y][gindex]
   else # if model is LA or RCG
       return vcat(pipie_cons(Y)[gindex],etavec(Y))
   end
end

# This function computes the test statistic
function kstesstat(ghat,G,Omegadiag,taun,solution)
    if sum(isnan.(ghat))>0
        return -100
    end
    dr,dg=size(G)
    KS=Model(Gurobi.Optimizer)
    #set_optimizer_attribute(KS,"outlev",0)
    @variable(KS,etavar[1:dg]>=taun/dg) #taun is a tuning parameter
    @objective(KS,Min,sum((sqrt(Omegadiag[r])*(ghat[r]-sum(G[r,l]*etavar[l] for l in 1:dg)))^2 for r in 1:dr))
    JuMP.optimize!(KS)
    if solution==true
        return G*value.(etavar)
    else
        return objective_value(KS)
    end
end

# This function computes the bootstrap statistic given the seed for a given frame
function ksbootseed(seed)
  Xb=X[rand(MersenneTwister(seed),1:N,N),:]
  ghatb=estimateg(Xb,gindex)
  return kstesstat(ghatb-ghat+etahat,G,Omegadiag,taun,false)
end

# This function computes the bootstrap statistic given the seed for the power simulations
function ksbootseedpower(seed, X, ghat, etahat, taun)
  Xb=X[rand(MersenneTwister(seed),1:size(X,1),size(X,1)),:]
  ghatb=estimateg(Xb,gindex)
  return kstesstat(ghatb-ghat+etahat,G,Omegadiag,taun,false)
end

# This function computes the bootstrap statistic given the seed for the model
# with stabel preferences
function ksbootseedstable(seed)
  rng1=MersenneTwister(seed)
  rng2=MersenneTwister(seed^2)
  rng3=MersenneTwister(2*seed+5)
  X1b=X1[rand(rng1,1:N1,N1),:]
  X2b=X2[rand(rng2,1:N2,N2),:]
  X3b=X3[rand(rng3,1:N3,N3),:]
  ghat1b=estimateg(X1b,gindex)
  ghat2b=estimateg(X2b,gindex)
  ghat3b=estimateg(X3b,gindex)
  ghatb=vcat(ghat1b,ghat2b,ghat3b)
  return kstesstat(ghatb-ghat+etahat,G,Omegadiag,taun,false)
end

function preferences(dYu,EU=false)
  U=collect(permutations(vec(1:dYu))) # All preference orders
  if EU
    Lot=[[] for i=1:dYu]
    Lot[1]=[0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5]
    Lot[2]=[0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    Lot[3]=[0.25, 0.0, 0.25, 0.0, 0.0, 0.25, 0.25]
    Lot[4]=[0.25, 0.2, 0.0, 0.15, 0.0, 0.0, 0.4]
    Lot[5]=[0.0, 0.2, 0.25, 0.15, 0.0, 0.25, 0.15]
    if dYu==6
        Lot[6]=[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    end
    St=zeros(length(U))
    A=zeros(length(Lot)-1,length(Lot[1]))
    B=[[] for i=1:length(U)]
    t=1
    for j in 1:length(U)
      M=Lot[U[j]]
      A=zeros(length(Lot)-1,length(Lot[1]))
      for k in 1:size(A,1)
        A[k,:]=M[k+1]-M[k]
      end
      KS=Model(Gurobi.Optimizer)
      set_optimizer_attribute(KS,"OutputFlag",0)
      dr,dg=size(A)
      @variable(KS,etavar[1:dg])
      @constraint(KS,cons[r=1:dr],sum(A[r,l]*etavar[l] for l in 1:dg)>=0.01 )
      @objective(KS,Min,etavar[1]^2); JuMP.optimize!(KS)
      if termination_status(KS)==MOI.OPTIMAL
        B[t]=U[j]
        t=t+1
      end
    end
    return B[[length(B[i])>0 for i in 1:length(B)]]
  else
    return U
  end
end
