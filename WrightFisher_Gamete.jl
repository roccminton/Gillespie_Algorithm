#=
Implementation of a simple Wright Fisher model with constant
population size

The population state has three variables
    "healthy" -> the number of healthy gametes
    "ill_ma" -> the number of gametes that are ill due to a incompatible mating
    "ill_mu" -> the number of gametes that are ill due to additional mutations
=#

using Distributions
using SparseArrays
using Random

healthy_pop(K) = Dict(
    "healthy" => K,
    "ill_mu" => 0,
    "ill_ma" => 0
    )

function WrightFisher.setupparameter(par, n0, historylength)
    return (
    par...,
    rndm = Vector{Int}(undef, 2),
    Mutation = Bernoulli(exp(-par.μ)),
    MutationLocation = 1:par.Nloci,
    traits = [inittraits(par, n0),inittraits(par, n0)],
    indices = Dict(
        "healthy" => [
            collect(1:n0["healthy"]),
            collect(1:n0["healthy"]),
        ],
        "ill_ma" => [
            collect(n0["healthy"]+1:n0["healthy"]+n0["ill_ma"]),
            collect(n0["healthy"]+1:n0["healthy"]+n0["ill_ma"])
        ],
        "ill_mu" => [
            collect(n0["healthy"]+n0["ill_ma"]+1:n0["healthy"]+n0["ill_ma"]+n0["ill_mu"]),
            collect(n0["healthy"]+n0["ill_ma"]+1:n0["healthy"]+n0["ill_ma"]+n0["ill_mu"])
        ],
    ),
    historylength = historylength,
    )
end

emptytraits(par) = spzeros(Bool,par.Nloci)

function inittraits(par, n0)
    #Generate healthy population
    traits = [emptytraits(par) for _ in 1:par.K]
    #if initial state is not perfectly healthy, one cannot accurately distribute the mutations
    if n0["ill_mu"] > 0 || n0["ill_ma"] > 0
        error("Cannot construct an initial population with the given population state")
    end
    return traits
end

function WrightFisher.birth!(ps, par, indx, old, new)
    #choose two gametes to mate
    rand!(par.rndm, par.indices["healthy"][old])
    #add a copy of the first parent to the population
    par.traits[new][indx] .= par.traits[old][par.rndm[1]]
    #check if they are compatible
    if iszero(sum(par.traits[old][par.rndm[1]] .* par.traits[old][par.rndm[2]]))
        #is a propagable mating
        #check for mutations at birth
        if rand(par.Mutation)
            #no mutation
            #and save the index in the approriate list
            ps["healthy"] += 1
            WrightFisher.saveorpush!(par.indices["healthy"][new], ps["healthy"], indx)
        else
            #mutation
            #chekc the location of the mutation
            mutation_location = rand(par.MutationLocation)
            #add the mutation to the gamete
            par.traits[new][indx][mutation_location] = true
            #check if one of the parental gametes already have a mutation at that position
            if par.traits[old][par.rndm[1]][mutation_location] || par.traits[old][par.rndm[2]][mutation_location]
                #mutation leads to incompatibility
                ps["ill_mu"] += 1
                WrightFisher.saveorpush!(par.indices["ill_mu"][new], ps["ill_mu"], indx)
            else
                #mutation does not lead to incompatibility
                ps["healthy"] += 1
                WrightFisher.saveorpush!(par.indices["healthy"][new], ps["healthy"], indx)
            end
        end
    else
        #is not a propagable mating
        ps["ill_ma"] += 1
        WrightFisher.saveorpush!(par.indices["ill_ma"][new], ps["ill_ma"], indx)
    end
end

const_popsize(popsize,par,t) = par.K
random_fluctuations(popsize,par,t) = round(Integer,rand(par.fluctuations))




#---
