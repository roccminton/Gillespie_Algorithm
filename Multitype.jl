"""
    module Multitype

    Implementation of a Logistic Model with a finite number of types.
    Mutation are possible at an independent rate μ and are specified by a
    transition matrix from one type to the next.
"""
module Multitype

export rates!, execute!

include("MainFunctions.jl")
import .Gillespie: chooseevent

include("BirthDeath.jl")
import .BirthDeath: LogisticRates!

function birth!(ps, i, pr)
    #Birth with or without mutation
    if rand() ≤ pr.μ
        #mutate to new type/species and add to species
        index_new_type = mutate(i,pr)
        #setup the size of the new type
        ps[index_new_type] += one(eltype(ps))
    else
        ps[i] += one(eltype(ps))
    end
    nothing
end

mutate(i, pr) = Gillespie.chooseevent(view(pr.M,:,i),1.0)

function death!(ps, i)
    ps[i] -= one(eltype(ps))
    nothing
end

function rates!(rates::Vector,ps,pr)
    @fastmath n = sum(ps)
    BirthDeath.LogisticRates!(rates,n,pr)
end

function rates!(rates::Matrix,ps,pr)
    rates[1,:] .= broadcast(x->BirthDeath.linearbirth(x,pr),ps)
    rates[2,:] .= broadcast(x->logisticdeath(x,ps,pr),ps)
end

logisticdeath(i,ps,pr) = ps[i]*(pr.death+ sum!(ps .* view(pr.competition,:,i)))

function execute!(i,ps,pr)
    if i==1
        birth!(ps,Gillespie.chooseevent(ps,sum(ps)),pr)
    elseif i==2
        death!(ps,Gillespie.chooseevent(ps,sum(ps)))
    else
        error("Index Error: No event #$i")
    end
end

end #end of module Multitype
