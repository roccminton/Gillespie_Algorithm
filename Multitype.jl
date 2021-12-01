#=
    module Multitype

    Implementation of a Logistic Model with a finite number of types.
    Mutation are possible at an independent rate μ and are specified by a
    transition matrix from one type to the next.

    Rates are represended as a Matrix, where the frist row gives the birh
    and the second row the death rate of each individual at the column index.
=#
module Multitype

export rates!, execute!, rungillespie

include("MainFunctions.jl")
import .Gillespie

using NamedTupleTools

setuprates(birthrate::Vector) = Matrix{eltype(birthrate)}(undef,(2,length(birthrate)))
setuphistory(time,n₀::Vector) = zeros(eltype(n₀),(length(time),length(n₀)))

function rungillespie(time,n₀,par;rescaled=false)
    #setup empty population history
    population_history = setuphistory(time,n₀)
    if rescaled
        Gillespie.run_gillespie!(
            time,n₀,
            generaterescaledparam(par),
            execute!,rates!,
            setuprates(par.birth),
            population_history
            )
    else
        Gillespie.run_gillespie!(
            time,n₀,
            (par...,diff=one(eltype(n₀))),
            execute!,rates!,
            setuprates(par.birth),
            population_history
            )
    end
    return population_history
end

function generaterescaledparam(par)
    par.birth .*= par.K
    par.death .*= par.K
    par.competition .*= par.K

    if haskey(par,:μ)
        fields = (field for field in fieldnames(par) if field != :μ)
        values = (par[field] for field in fields)
        return NamedTuple{(fields...,:μ,:diff)}(values...,par.μ/par.K,inv(par.K))
    else
        return (par...,diff=inv(par.K))
    end
end

function birth!(ps, i, pr)
    #Birth with or without mutation
    if rand() ≤ pr.μ
        #mutate to new type/species and add to species
        index_new_type = mutate(i,pr)
        #setup the size of the new type
        ps[index_new_type] += pr.diff
    else
        ps[i] += pr.diff
    end
    nothing
end

mutate(i, pr) = Gillespie.chooseevent(view(pr.M,:,i),1.0)

function death!(ps, i,pr)
    ps[i] -= pr.diff
    nothing
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
        death!(ps,div(i,2),pr)
    end
end

end #end of module Multitype
