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

#=
Rates function if birth and death rates are equal for all individuals.
In that case the rates Vector has only two elements.
=#
function rates!(rates::Vector,ps,pr)
    @fastmath n = sum(ps)
    BirthDeath.LogisticRates!(rates,n,pr)
end

#=
Rates function if birth and death rates vary from one individual to another.
In that case the rates are represented as a Matrix. The first row holds the
birth rates the second row holds the death rates. Individuals are listed colomn-
wise.
=#

function rates!(rates::Matrix,ps,pr)
    #linear birth
    view(rates,1,:) .= ps .* pr.birth
    #for each individual get a total death rate of
    #   n_x (d + Σ c(x,y)n_y)
    #death rate
    view(rates,2,:) .= pr.death
    #additional death through competition
    view(rates,2,:) .+= pr.competition * ps
    #total death rate
    view(rates,2,:) .*= ps
end

function execute!(i,ps,pr)
    if isodd(i)
        birth!(ps,div(i+1,2),pr)
    else
        death!(ps,div(i,2))
    else
        error("Index Error: No event #$i")
    end
end

end #end of module Multitype
