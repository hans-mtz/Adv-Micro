# Pkg.add("JuMP")
# using JuMP
# # On Mac, this might be
# ENV["GUROBI_HOME"] = "/Library/gurobi911/mac64"
# ENV["GRB_LICENSE_FILE"]="/Library/gurobi911/mac64/gurobi.lic"
# # import Pkg
# Pkg.add("Gurobi")
# Pkg.build("Gurobi")
# Pkg.add("MathProgBase")
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("DataFrames")

using LinearAlgebra

using MathProgBase
using JuMP, Gurobi
#using Clp
using DataFrames
using CSV

## Setting-up directory
# tempdir1=@__DIR__
# repdir=tempdir1[1:findfirst("ReplicationAK",tempdir1)[end]]
# appname="FirstApp"
# rootdir=repdir*"/"*appname
# diroutput=repdir*"/Output_all"
# dirdata=repdir*"/Data_all"

## Setting directory Hans
# cd("/Volumes/SSD Hans/Github/TopicsDecisionMaking/ReplicationAK")
rootdir="/Volumes/SSD Hans/Github/Adv Micro/Adv-Micro/PS1"
dirdata="/Volumes/SSD Hans/Github/Adv Micro/Adv-Micro/PS1/Data_all"


## Functions
include(rootdir*"/ED_det_test.jl")  # ED deterministic test function
include(rootdir*"/ED_data_load.jl") # Function that loads the data

## Testing singles
rho,cve=ED_data_load(dirdata,"singles") # Data loading
rate_singles=ED_det_test(rho,cve,stepdum) # Testing

## Loading and parameters
#n=sample size, T=time periods, K= goods
n,T,K=size(rho)
 # Data loading
rho,cve=ED_data_load(dirdata,"couples")
stepdum= .05

## testing
include(rootdir*"/ED_det_test-hans.jl")
rate_couples, test=ED_det_test(rho,cve,stepdum)
