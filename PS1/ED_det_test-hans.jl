function ED_det_test(rho,cve,stepdum)
    n,T,K=size(rho)
    dummy=Model(Gurobi.Optimizer)
    set_optimizer_attribute(dummy,"OutputFlag",0)
    @variable(dummy,v[1:T])
    @constraint(dummy,ineq,v[1]>=0)
    @objective(dummy,Min,1)
    optimize!(dummy)
    ##no solution type
    val=primal_status(dummy)


    deltat=collect(.1:stepdum:1)
    soldet=zeros(size(deltat,1),size(deltat,1),n)
## Model initialization
    ind=1
    p=rho[ind,:,:]'
    q=cve[ind,:,:]'
    δA=.9
    δB=.9
    # n=3

    EDtest=Model(Gurobi.Optimizer)
    set_optimizer_attribute(EDtest,"OutputFlag",0)
    @variable(EDtest,vA[1:T])
    @variable(EDtest,vB[1:T])
    @constraint(EDtest,ineq[t=1:T,s=1:T],δA^(t-1)*(vA[t] - vA[s]) + δB^(t-1)*(vB[t] - vB[s])>=p[:,t]'*(q[:,t]-q[:,s]))
## loop

    for ind=1:n
        p=rho[ind,:,:]'
        q=cve[ind,:,:]'
        for dd=1:size(deltat,1),db=1:size(deltat,1)
            # for i=1:T, j=1:T
            #     delete(EDtest,ineq[i,j])
            # end
            δA=deltat[dd]; δB=deltat[db]
            # @constraint(EDtest,ineq[t=1:T,s=1:T],δA*(vA[t] - vA[s]) + δB*(vB[t] - vB[s])>=p[:,t]'*(q[:,t]-q[:,s]))
            for t=1:T, s=1:T
                for a=1:T
                    for b=1:T
                        if a==b
                            continue
                        end
    ##
                          # a=1;b=3;δA=0.8;δB=0.7
                          set_normalized_coefficient(ineq[a,b],vA[a],δA)
                          set_normalized_coefficient(ineq[a,b],vA[b],-δA)
                          set_normalized_coefficient(ineq[a,b],vB[a],δB)
                          set_normalized_coefficient(ineq[a,b],vB[b],-δB)
                          # println(ineq[a,b])
                      end
                  end
 # #

                  set_normalized_rhs(ineq[t,s], p[:,t]'*(q[:,t]-q[:,s]))
                  # println("Household ",ind,". ",ineq[t,s])
            end


            @objective(EDtest,Min,1)
            optimize!(EDtest)
            soldet[dd,db,ind]=primal_status(EDtest)==val
        end
        ind % 100 == 0 && println("Household ",ind)
    end

    ## Rejection Rate
    rate=1-sum((sum(soldet,dims=(1,2))[:].>=ones(n))*1)/n
return rate, soldet
end


# list_of_constraint_types(EDtest)
