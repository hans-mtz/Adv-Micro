using Distributed
using Statistics
using DataFrames, CSV
addprocs(7)

@everywhere begin
  using Random
  using Combinatorics
  using LinearAlgebra
  using JuMP
  using Gurobi
  using KNITRO
end
@everywhere model=$(ARGS[1])   # put "LA", "RCG", or "RUM"
datatype=ARGS[2]               # put "high","low", or "medium"
ExpUt=("EU"=="EU")             # put true or false

# @everywhere model="RUM"   # put "LA", "RCG", or "RUM"
# datatype="pooled"               # put "pooled", "high","low", or "medium"/frame
# ExpUt=("EU"=="EU")          # put true or false

println(model)
println(datatype)
println(ExpUt)

## Defining the file directories
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("ConsiderationSetExperiment",tempdir1)[end]]
dirresults=rootdir*"/results/"

## Functions
@everywhere include($(rootdir)*"/Codes/Julia_code_nail/newfunctions/functions_common_new.jl")

## Common parameters
dYm=5                               # Number of varying options in menus
model=="RUM" ? dYu=6 : dYu=5        # For RUM there are 6 options instead of 5
Menus=collect(powerset(vec(1:dYm))) # Menus
gindex=gindexes(Menus,dYu)          # Indices that correspond to nonzero linearly independent frequencies
U=preferences(dYu,ExpUt)            # All preference orders
G=matrixcons(gindex, Menus, U)      # Marix of 0 and 1
Omegadiag=ones(size(G,1))           # Vector of weigts

## Data
X=Matrix(CSV.read(rootdir*"/Codes/Julia_code_nail/newfunctions/data/menu_choice_$(datatype).csv", DataFrame))
# Sample size
N=size(X,1)
# Smallest sample per menu
Ntau=nintaun(X)
# Tuning parameter
taun=sqrt(log(Ntau)/Ntau)

## Testing
println("Testing...")
# Estimates of g
ghat=estimateg(X,gindex)
etahat=kstesstat(ghat,G,Omegadiag,taun,true)
@everywhere begin
  X=$X; N=$N; gindex=$gindex; ghat=$ghat; G=$G
  Omegadiag=$Omegadiag; taun=$taun; etahat=$etahat; Menus=$Menus
end
# Test statistic
Tn=N*kstesstat(ghat,G,Omegadiag,0.0,false)
# Bootstrap statistics
Boot=pmap(ksbootseed,1:1000) # use ksbootseed function
# Pvalue. If it is 0.0, then pvalue<0.001
pvalue=mean(Tn.<=N.*collect(Boot))

## Saving the output
if ExpUt
  CSV.write(rootdir*"/Codes/Julia_code_nail/newfunctions/results/Tn_new_pval_EU-$(model)_$(datatype).csv", DataFrame(Tn=Tn, pvalue=pvalue))
else
  CSV.write(rootdir*"/Codes/Julia_code_nail/newfunctions/results/Tn_new_pval_$(model)_$(datatype).csv", DataFrame(Tn=Tn, pvalue=pvalue))
end
println("Done!")
