using Distributed
using Statistics
using DataFrames, CSV
using Plots
addprocs(4)
# Pkg.add("Combinatorics")

@everywhere begin
  using Random
  using Combinatorics
  using LinearAlgebra
  using JuMP
  using Gurobi
  ##using KNITRO
  using Ipopt
end

## Defining the file directories
# tempdir1=@__DIR__
# rootdir=tempdir1[1:findfirst("ABKKChoiceExperiment",tempdir1)[end]]
rootdir="/Volumes/SSD Hans/Github/Adv Micro/Adv-Micro/Midterm/ABKKChoiceExperiment"
cd(rootdir)
# Functions
@everywhere include($(rootdir)*"/newfunctions/functions_common_new.jl")

## Data for all 3 frames
X1=Matrix(CSV.read(rootdir*"/data/menu_choice_high.csv", DataFrame))
X2=Matrix(CSV.read(rootdir*"/data/menu_choice_medium.csv", DataFrame))
X3=Matrix(CSV.read(rootdir*"/data/menu_choice_low.csv", DataFrame))
println("Data is ready!")
# Sample sizes
N1=size(X1,1)
N2=size(X2,1)
N3=size(X3,1)
# Smallest sample size
N=minimum([N1,N2,N3])
# Smallest sample per menu
Ntau=nintaun(X1)
# Tuning parameter as suggested in KS
taun=sqrt(log(Ntau)/Ntau)


## Re-Run when changing models
@everywhere model="LA"   # put "LA", "RCG", or "RUM"
println(model)
ExpUt=true          # put true if expected utility  or false if general
# println(ExpUt)
# Common parameters
dYm=5                               # Number of varying options in menus
model=="RUM" ? dYu=6 : dYu=5        # For RUM there are 6 options instead of 5
Menus=collect(powerset(vec(1:dYm))) # Menus
gindex=gindexes(Menus,dYu)          # Indices that correspond to nonzero linearly independent frequencies

## Victor's example
U=preferences(dYu, ExpUt)           # All preference orders
G=matrixcons(gindex, Menus, U)      # Matrix of 0 and 1
G=vcat(G,G,G)                       # To impose stability we need to repeat G for every frame
Omegadiag=ones(size(G,1))           # Vector of weigts


## Testing
println("Testing...")
# Estimates of g for all 3 frames
ghat1=estimateg(X1,gindex)
ghat2=estimateg(X2,gindex)
ghat3=estimateg(X3,gindex)
ghat=vcat(ghat1,ghat2,ghat3)
etahat=kstesstat(ghat,G,Omegadiag,taun,true)
@everywhere begin
  X1=$X1; N1=$N1; X2=$X2; N2=$N2; X3=$X3; N3=$N3; N=$N
  gindex=$gindex; ghat=$ghat; G=$G;
  Omegadiag=$Omegadiag; taun=$taun; etahat=$etahat; Menus=$Menus
end
# Test statistic
Tn=N*kstesstat(ghat,G,Omegadiag,0.0,false)
# Bootstrap statistics
Boot=pmap(ksbootseedstable,1:1000) # use ksbootseedstable function
# Pvalue. If it is 0.0, then pvalue<0.001
pvalue=mean(Tn.<=N.*collect(Boot))

## optional
##histogram(N.*collect(Boot), bins = :scott, weights = repeat(1:5, outer = 200))


## Saving the output
if ExpUt
  CSV.write(rootdir*"/results/Tn_new_pval_EU-$(model)_stable.csv", DataFrame(Tn=Tn, pvalue=pvalue))
else
  CSV.write(rootdir*"/results/Tn_new_pval_$(model)_stable.csv", DataFrame(Tn=Tn, pvalue=pvalue))
end
println("Done!")
