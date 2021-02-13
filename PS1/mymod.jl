# module mymod

## Getting mean ρ⋆c⋆ for each individual
    function avg_inc(ind::Integer,rho,cve )
        n,T,K=size(rho)
        inc=0
        # ind=1
        for t=1:T
            # t=4
            p=rho[ind,t,:]'
            q=cve[ind,t,:]'

            pq=p*q'
            # println("income for individual ", ind," in period ",t," is ",pq)
            inc=inc+pq
        end
        inc=inc/4
        # print("average income is ",inc)
        w=[inc-rho[ind,t,:]'*cve[ind,t,:] for t=1:4]

        return inc, w'
        # p=rho[ind,t,:]'
        # q=cve[ind,t,:]'
        # w[t]=inc-p*q'
    end

    function avg_inc2d(rho,cve )
        T,K=size(rho)
        inc=0
        # ind=1
        for t=1:T
            # t=4
            p=rho[t,:]'
            q=cve[t,:]'

            pq=p*q'
            # println("income for simulation ", t," is ",pq)
            inc=inc+pq
        end
        inc=inc/T
        # print("average income is ",inc,"\n")
        w=[inc-rho[t,:]'*cve[t,:] for t=1:T]

        return inc, w'
        # p=rho[ind,t,:]'
        # q=cve[ind,t,:]'
        # w[t]=inc-p*q'
    end

    function avg_inc_v(rho,cve )
        K,T=size(rho)
        # inc=0
        # # ind=1
        # for t=1:T
        #     # t=4
        #     p=rho[t,:]'
        #     q=cve[t,:]'
        #
        #     pq=p*q'
        #     # println("income for simulation ", t," is ",pq)
        #     inc=inc+pq
        # end
        # inc=inc/T
        # print("average income is ",inc,"\n")
        inc=[rho[:,t]'*cve[:,t] for t=1:T]

        return inc
        # p=rho[ind,t,:]'
        # q=cve[ind,t,:]'
        # w[t]=inc-p*q'
    end
    ## Simulations for one individual all time periods

    function sim(ind,rho,cve,h,tol,γ)
        n,T,K=size(rho)
        # ind=2
        # t=c_1
        # h=0
        # tol=1e-3

        p=rho[ind,:,:]'
        q=cve[ind,:,:]'
        inc, w = avg_inc(ind,rho,cve)
        σ=maximum(abs.(w))

        # ub=inc+σ
        # lb=inc-σ

        # uσ=maximum(w)
        # lσ=minimum(w)
        S=T

        Random.seed!(h*13+17)
        s_s=rand(S,K)

        # s=SobolSeq(S)
        # skip(s, K*h*10)
        # s_s=hcat([next!(s) for i = 1*h*10:K*h*10]...)
        # w_s=[lσ+s_s[i,j]*(uσ-lσ) for i=1:S, j=1:K]
        # w_r=[-1+s_s[i,j]*2 for i=1:S, j=1:K]
        # c_t=[lb+s_s[i,j]*(ub-lb) for i=1:S, j=1:K]

        # p_s=repeat(p,S)

        # w_t=[-σ+s_s[i,j]*2*σ for i=1:S, j=1:K]'
        w_t=[-1*γ+s_s[i,j]*2*γ for i=1:S, j=1:K]'

        c_1=q-w_t
        # c_2=repeat(q,S)-w_s
        # c_3=repeat(q,S)-w_r

        inc_h=avg_inc_v(p,c_1)
        inc_t=avg_inc_v(p,q)
        # avg_inc2d(p_s,c_t)
        # avg_inc2d(p_s,c_2)
        # avg_inc2d(p_s,c_3)
        iah=sum(inc_h)/4
        inc_t-inc_h.<=tol
        # iah,m=avg_inc2d(p,c_1)
        # iat,n=avg_inc2d(p,q)
        # iat==inc

        A=zeros(T,T)
        B=ones(T,T)
        C=[(p[:,i]-p[:,j])'*(c_1[:,i]-c_1[:,j]) for i=1:T,j=1:T]
        vals,vecs=eigen(C)

        # false || println("income not OK ",iat," got ",iah," Diff=",iat-iah); return 0
        #
        # true && println("Hard ineq OK!"); return 1

        v=zeros(T)

        # true & true && println("NSD OK!"); return 1
        # if isreal(vals)
        #     cum=0
        #     for i=1:size(vals)
        #         eigs=vals[i]<= 0
        #         cum=cum+eigs
        #     end
        #     cum>=4 && return 1
        # end
        println("Evaluating individual ",ind," sim=",h)
        # if abs(inc-iah)<=tol
        #     println("income not OK ",inc," got ",iah," Diff=",iat-iah)
        #     return 0
        if sum(A.>=C)==sum(B)
            println("Hard ineq OK!")
            return 1
        elseif isreal(vals) & sum(vals.<v)>=4
            println("NSD OK!")
            return 1
        else
            println("Nothing for you, bro, sorry");
            return 0
        end
    end
##
# end  # module mymod
