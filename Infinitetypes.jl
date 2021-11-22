#=
    module InfiniteTypes

    Implementation of a Logistic Model with a infinite number of possible types.
    Mutation are possible at an independent rate μ and are specified by a
    transition kernel from one type to the next.

    Rates are represended as a Dictionary of Vectors, where the frist entry
    of the Vector gives the birh and the second entry the death rate of the
    individual type with that key.
=#

module InfiniteTypes

include("Onetype.jl")
import .OneType


export rates!, execute!

function rates!(rates::Dict,ps::Dict,pr)
    for (trait,rate) in ps
        nₓ = ps[trait][1]
        #birthrate n_x * b(x)
        rate[1] = nₓ*ps[trait][2]
        #deathrate n_x * (d(x) + Σ c(x,y) n_y)
        rate[2] = nₓ* ps[trait][3]
        for (traittuple,c) in pr.compdict
            t₁,t₂ = traittuple
            t₁ == trait && rate[2] += nₓ * ps[t₂][1] * c
        end
    end
end

function birth!(ps, rates, pr, trait)
    #Birth with or without mutation
    if rand() ≤ pr.μ
        #mutate to new type/species and add to species
        new_trait = mutate(trait,pr)
        #setup the size of the new type
        if haskey(ps,new_trait)
            ps[new_trait][1] += pr.diff
        else
            addnewtrait!(ps,rates,pr)
        end
    else
        ps[trait] += pr.diff
    end
    nothing
end

function death!(ps,trait,pr)
    ps[trait] -= pr.diff
end

function execute!()

function addnewtrait!(ps,rates,pr, trait)
    #add to population state
    ps[trait] = [pr.diff,pr.birth(trait),pr.death(trait)]
    #initialize to rates
    rates[trait] = valtype(rates)(undef,2)
    #set competition
    for other_trait in keys(ps)
        pr.compdict[(trait,other_trait)] = pr.competition(trait,other_trait)
        pr.compdict[(other_trait,trait)] = pr.competition(other_trait,trait)
    end
end

end  # end of module InfiniteTypes
