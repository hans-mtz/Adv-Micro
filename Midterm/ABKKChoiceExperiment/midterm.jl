using LinearAlgebra
using Plots
using TexTables
using Latexify



# σ=0.5
## Utility function
# Takes in values of sigma(σ), utility dimensions(number of loteries=dyu) and in the case of RUM considers
# to add the default option in (m=true)
# This function returns the raking for a given σ.
function crra(σ, m, dyu)
    Lot=[[] for i=1:dyu]
    Lot[1]=[0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5]
    Lot[2]=[0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0]
    Lot[3]=[0.25, 0.0, 0.25, 0.0, 0.0, 0.25, 0.25]
    Lot[4]=[0.25, 0.2, 0.0, 0.15, 0.0, 0.0, 0.4]
    Lot[5]=[0.0, 0.2, 0.25, 0.15, 0.0, 0.25, 0.15]
    if dyu==6
        Lot[6]=[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    end
    prize=[50, 48, 30, 14, 12, 10, 0]
    if m==true
        ix=[true, true, true, true, false, true, true]
        Lot=[Lot[i][ix] for i=1:length(Lot)]
        prize=prize[ix]
    end
    if σ==1
        u=log.(prize[1:end-1])
    else
        u=(prize[1:end-1].^(1-σ))./(1-σ)
    end
    eu=[dot(vec(Lot[i][1:end-1]),u) for i=1:dyu ]
    r=sort(eu,rev=true)
    i=indexin(eu,r)
    return i
end

##
# unique(x -> crra(x), -1:0.1:1)
#
# ranking=crra(0.5)
# Menus[31]
# indexin(Menus[27],ranking)
##

# function get_B(M,U)
#     sigma=[i for i=-1:0.05:1]
#     U=(crra.(sigma))
#     U=unique(U)
#     B=zeros(length(gind),length(U))
#     for j=1:length(M)
#         if length(M[j])>0
#             for i=1:length(U)
#                 if i==minimum(collect(indexin(M[j],U)))
#                     B[j,i]=1.0
#                 # B[j,i]=i==minimum(collect(indexin(M[j],ranking)))
#         end
#     end
#     return B
# end
#
# ##
#
# B=get_B(Menus,σ)

##

# rs=[crra(i) for i=-1:0.01:1]
#
# unique(collect(rs))

## Getting rankings for σ ∈ [-1,1] in steps of 0.01

sigma=[i for i=-1:0.01:1]
Uc=(crra.(sigma,false,dYu))
Ucn=unique(Uc) #Identifying the unique different rankings (only 6 of them for CRRA)
s_grid=sigma[indexin(Ucn,Uc)] #Saving the different values of σ that provide unique rankings

##
Gc=matrixcons(gindex, Menus, Ucn)   # Using VK's code to get G matrix
G=vcat(Gc,Gc,Gc)                    # To impose stability we need to repeat G for every frame
Omegadiag=ones(size(G,1))           # Vector of weigts


##

# Computing G
# U is the set of preferences
# M is the set of menus
# gindexsh are coordinates of nonzero linearly independent p_pi
# function matrixconsCRRA(gindexsh, M, U)
#     dYu=length(U[1])
#     dYm=length(M)
#     if dYu==6
#         M=[vcat(1,M[i].+1) for i in eachindex(M)]
#         M[1]=[]
#     end
#     d=length(U)
#     d2=length(gindexsh)
#     B=zeros(d2,d)
#     m1=1
#     for j in 1:d
#         pp=zeros(dYm,dYu)
#         for i in eachindex(M)
#             if length(M[i])>0 # Skipping empty menu
#                 for k in 1:dYu
#                     if k==M[i][argmax(U[j][M[i]])]
#                         pp[i,k]=1.0
#                     end
#                 end
#             end
#         end
#         # println(pp)
#         B[:,m1]=pp[gindexsh] #picking only relevant elements, indexes that are not always zero
#         m1=m1+1
#     end
#     if dYu==6 # If model is RUM then return B
#         return B
#     else #otherwise return
#         return [[B zeros(size(B,1),dYm)];[zeros(dYm,size(B,2),) I]]
#     end
# end


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
Tn_CRRA_RUM=N*kstesstat(ghat,G,Omegadiag,0.0,false)
# Bootstrap statistics
Boot=pmap(ksbootseedstable,1:100) # use ksbootseedstable function
# Pvalue. If it is 0.0, then pvalue<0.001
pvalue_CRRA_RUM=mean(Tn_CRRA_RUM.<=N.*collect(Boot))

# Making histogram and saving it
histogram(collect(N.*Boot))#, bins = :scott, weights = repeat(1:5, outer = 200))
title!("CRRA, RUM")
savefig("/Volumes/SSD Hans/Github/Adv Micro/Adv-Micro/Midterm/ABKKChoiceExperiment/figs/CRRA_RUM")

## Getting rankings for σ ∈ [-1,1] in steps of 0.01
#LA model
sigma=[i for i=-1:0.01:1]
Uc=(crra.(sigma,false,5))
Ucn=unique(Uc) #Identifying the unique different rankings (only 6 of them for CRRA)
s_grid=sigma[indexin(Ucn,Uc)] #Saving the different values of σ that provide unique rankings

##
Gc=matrixcons(gindex, Menus, Ucn)   # Using VK's code to get G matrix
G=vcat(Gc,Gc,Gc)                    # To impose stability we need to repeat G for every frame
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
Tn_CRRA_LA=N*kstesstat(ghat,G,Omegadiag,0.0,false)
# Bootstrap statistics
Boot=pmap(ksbootseedstable,1:100) # use ksbootseedstable function
# Pvalue. If it is 0.0, then pvalue<0.001
pvalue_CRRA_LA=mean(Tn_CRRA_LA.<=N.*collect(Boot))

## Making histogram and saving it
histogram(collect(N.*Boot))#, bins = :scott, weights = repeat(1:5, outer = 200))
title!("CRRA, LA")
savefig(rootdir*"/figs/CRRA_LA")

## Saving the output

results=DataFrame(  Tn     =[Tn_CRRA_LA,Tn_CRRA_RUM],
                    pvalue =[pvalue_CRRA_LA,pvalue_CRRA_RUM])
CSV.write("./results/results.csv",
            DataFrame(Tn=[Tn_CRRA_LA,Tn_CRRA_RUM],
            pvalue=[pvalue_CRRA_LA,pvalue_CRRA_RUM]))
println("Done!")

mdtable(results)

## Question 2

# get distributions under LA using CRRA, using High and Medium cost treatments

## Getting rankings for σ ∈ [-1,1] in steps of 0.01
#LA model
sigma=[i for i=-1:0.01:1]
Uc=(crra.(sigma,false,5))
Ucn=unique(Uc) #Identifying the unique different rankings (only 6 of them for CRRA)
s_grid=sigma[indexin(Ucn,Uc)] #Saving the different values of σ that provide unique rankings

##
Gc=matrixcons(gindex, Menus, Ucn)   # Using VK's code to get G matrix
G=vcat(Gc,Gc)                    # To impose stability we need to repeat G for every frame
Omegadiag=ones(size(G,1))           # Vector of weigts


## Testing
println("Testing...")
# Estimates of g for all 3 frames
ghat1=estimateg(X1,gindex)
ghat2=estimateg(X2,gindex)
ghat3=estimateg(X3,gindex)
ghat=vcat(ghat1,ghat2)
etahat=kstesstat(ghat,G,Omegadiag,taun,true)


## Plots
plot(1:6,(G'*etahat)[1:6]./sum((G'*etahat)[1:6]),
            title="Low-Cost Risk Aversion Distribution",
            leg=:bottomright,
            label="Prediction CRRA-LA", line=(:dot,4),
            xlabel="σ ∈ [-1,1]") #Prediction
plot!(1:6,(Gc'ghat3)[1:6]./sum((Gc'*ghat3)[1:6]),label="Data", line=4)
xticks!(1:6,string.(s_grid))
savefig("../figs/pred_LA")  #Distribution for High and Medium

plot(1:6,(Gc'ghat1)[1:6]./sum((Gc'*ghat1)[1:6]), leg=:bottomright,
        label="High", title="Risk aversion distributions CRRA-LA", line=4,
        xlabel="σ ∈ [-1,1]") #Distribution for High and Medium
plot!(1:6,(Gc'ghat2)[1:6]./sum((Gc'*ghat2)[1:6]), label="Medium", line=4)
plot!(1:6,(Gc'ghat3)[1:6]./sum((Gc'*ghat3)[1:6]), label="Low", line=4)
xticks!(1:6,string.(s_grid))
savefig("../figs/dist_LA")
## Error LA prediction
low_pred_CRRA_LA=(G'*etahat)[1:6]./sum((G'*etahat)[1:6])
low_data=(Gc'ghat3)[1:6]./sum((Gc'*ghat3)[1:6])
err_LA_low=sum((low_data.-low_pred_CRRA_LA).^2)
err_LA_low=sqrt(err_LA_low)

# get distributions under RUM using CRRA, using High and Medium cost treatments

## Getting rankings for σ ∈ [-1,1] in steps of 0.01
#RUM model
sigma=[i for i=-1:0.01:1]
Uc=(crra.(sigma,false,dYu))
Ucn=unique(Uc) #Identifying the unique different rankings (only 6 of them for CRRA)
s_grid_RUM=sigma[indexin(Ucn,Uc)] #Saving the different values of σ that provide unique rankings

##
Gc=matrixcons(gindex, Menus, Ucn)   # Using VK's code to get G matrix
G=vcat(Gc,Gc)                    # To impose stability we need to repeat G for every frame
Omegadiag=ones(size(G,1))           # Vector of weigts


## Testing
println("Testing...")
# Estimates of g for all 3 frames
ghat1=estimateg(X1,gindex)
ghat2=estimateg(X2,gindex)
ghat3=estimateg(X3,gindex)
ghat=vcat(ghat1,ghat2)
etahat=kstesstat(ghat,G,Omegadiag,taun,true)

## Plots

plot((1:10,["aqui"]),(G'*etahat)[1:10]./sum((G'*etahat)[1:10]),
            title="Low-cost Riks Aversion Distribution",
            leg=:bottomright,
            xlabel="σ ∈ [-1,1]", label="Prediction CRRA-RUM", line=(:dot,4))
xticks!(1:10,string.(s_grid_RUM)) #Prediction
plot!(1:10,(Gc'ghat3)[1:10]./sum((Gc'*ghat3)[1:10]),label="Data", line=4)
savefig("./figs/pred_RUM") #Distribution for High and Medium
plot(1:10,(Gc'ghat1)[1:10]./sum((Gc'*ghat1)[1:10]), leg=:bottomright,
        label="High", title="Risk aversion distributions RUM-CRRA", line=4,
        xlabel="σ ∈ [-1,1]") #Distribution for High and Medium
plot!(1:10,(Gc'ghat2)[1:10]./sum((Gc'*ghat2)[1:10]), label="Medium", line=4)
plot!(1:10,(Gc'ghat3)[1:10]./sum((Gc'*ghat3)[1:10]), label="Low", line=4)
xticks!(1:10,string.(s_grid_RUM))
savefig("./figs/dist_RUM")
## Error LA prediction
low_pred_CRRA_RUM=(G'*etahat)[1:10]./sum((G'*etahat)[1:10])
low_data=(Gc'ghat3)[1:10]./sum((Gc'*ghat3)[1:10])
err_RUM_low=sum((low_data.-low_pred_CRRA_RUM).^2)
err_RUM_low=sqrt(err_RUM_low)

## save results

err_t=DataFrame(errors=[err_LA_low,err_RUM_low])
CSV.write("../results/pred_errs.csv",err_t)
