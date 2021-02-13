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
    soldet=zeros(size(deltat,1),n)
    ## Model initialization
    ind=1
    p=rho[ind,:,:]'
    q=cve[ind,:,:]'
    delta=.9
    EDtest=Model(Gurobi.Optimizer)
    set_optimizer_attribute(EDtest,"OutputFlag",0)
    @variable(EDtest,v[1:T])
    @constraint(EDtest,ineq[t=1:T,s=1:T],v[t]-v[s]>=delta^(-t+1)*p[:,t]'*(q[:,t]-q[:,s]))
    for ind=1:n
      for dd=1:size(deltat,1)
        delta=deltat[dd]
        p=rho[ind,:,:]'
        q=cve[ind,:,:]'
        for t=1:T
          for s=1:T
            set_normalized_rhs(ineq[t,s], delta^(-t+1)*p[:,t]'*(q[:,t]-q[:,s]))
          end
        end
        @objective(EDtest,Min,1)
        optimize!(EDtest)

        soldet[dd,ind]=primal_status(EDtest)==val
      end
    end

    ## Rejection Rate
    rate=1-sum((sum(soldet,dims=1)[:].>=ones(n))*1)/n
return rate
end
