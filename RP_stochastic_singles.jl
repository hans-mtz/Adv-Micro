## Problem 3

# Probability of individual i is rational P(i is rational=
# P(0≥ρ[t]'(c[t]-w[t]+c[s]+w[s]), E[ρ[t]'c[t]]-ρ⋆c⋆)=
# N^-1∑1{0≥ρ[t]'(c[t]-w[t]+c[s]+w[s]),E[ρ[t]'c[t]]-ρ⋆c⋆)}

## Packages
using LinearAlgebra
using MathProgBase
using JuMP, Gurobi
#using Clp
using DataFrames
using CSV
# Pkg.add("Sobol")
using Sobol
using Random
# using mymod




## Setting directory Hans and getting Data_all
# cd("/Volumes/SSD Hans/Github/TopicsDecisionMaking/ReplicationAK")
difun="/Volumes/SSD Hans/Github/TopicsDecisionMaking/ReplicationAK/FirstApp"
dirdata="/Volumes/SSD Hans/Github/TopicsDecisionMaking/ReplicationAK/Data_all"

include(rootdir*"/ED_data_load.jl") # Function that loads the data
rho,cve=ED_data_load(dirdata,"singles")
n,T,K=size(rho)

## Simulating for ind 2
include("mymod.jl")


tol=1e-1
rej_rat=0
Sims=1000
γ=3

for i=1:n
    i2=[sim(i,rho,cve,h,tol,γ) for h=1:Sims]
    rate=1-(sum(i2)/Sims)
    rej_rat=rej_rat+rate
end
rej_rat=rej_rat/n
println(rej_rat)

# i2=[sim(185,rho,cve,h,tol) for h=1:Sims]
# rate=1-(sum(i2)/Sims)

## deterministic
include("ED_det_WL.jl")
rate_couples=ED_det_test(rho,cve)
