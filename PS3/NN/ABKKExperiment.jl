# using Distributed
# using Statistics
using DataFrames, CSV
# using Plots
# addprocs(7)

# @everywhere begin
#   using Random
#   using Combinatorics
#   using LinearAlgebra
#   using JuMP
#   using Gurobi
#   ##using KNITRO
#   using
# end

# tempdir1=@__DIR__
tempdir1="/Volumes/SSD Hans/Github/Adv Micro/Adv-Micro/PS3/NN"
# rootdir=tempdir1[1:findfirst("TopicsDecisionMaking",tempdir1)[end]]
# rootdir=tempdir
# data=CSV.read(tempdir1*"/data/ABKK_full_experiment.csv", DataFrame)

dnn=CSV.read(tempdir1*"/dnn.csv", DataFrame)
