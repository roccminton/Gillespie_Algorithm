"""
Implementation of a simple Wright Fisher model with constant
population size
"""

module WrightFisher

using Random
using Distributions
using ProgressMeter
using SparseArrays

include("DiploidModel2.jl")
include("Plotting.jl")
include("MainFunctions.jl")

import .DiploidModel2
import .PlotFromDicts
import .Gillespie

#execute simulation
function run_wrightfisher(time, n₀, model_parameter)

    #setup empty population history
    l = length(time)
    population_history = Dict(x => zeros(valtype(n₀), l) for x in keys(n₀))

    #execute simulation
    mainiteration!(
        n₀,
        time,
        setupparameter(model_parameter, n₀, l),
        population_history,
    )

    return population_history
end

setupparameter(par, n0, historylength) = (
    par...,
    rndm = Vector{Int}(undef, 2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = DiscreteUniform(1, par.Nloci),
    traits = [
        DiploidModel2.inittraitsfromdict(par, n0),
        DiploidModel2.inittraitsfromdict(par, n0)
        ],
    indices = Dict(
        "healthy" => [
            collect(n0["Ill"]+1:n0["PopSize"]),
            collect(n0["Ill"]+1:n0["PopSize"]),
        ],
        "ill" => [
            collect(1:n0["Ill"]),
            collect(n0["Ill"]+1:n0["PopSize"])
        ],
    ),
    historylength = historylength,
)

function mainiteration!(n0, time, par, pop_hist)
    #run simulation
    @showprogress for t in time
        #save one step evolution
        Gillespie.saveonestep!(pop_hist, t, n0, par)
        #empty population state
        for k in keys(n0)
            n0[k] = 0
        end
        #execute one step of the simulation
        onestep!(n0, par, t)
    end
end

oldnew(t) = isodd(t) ? (1, 2) : (2, 1)

function onestep!(n0, par, t)
    old, new = oldnew(t)
    #create next generation
    for child = 1:par.popsize(n0,par,t)
        birth!(n0, par, child, old, new)
    end
    #cut the index lists to smaller size if necessary
    cuttosize!(par.indices["healthy"][new], n0["PopSize"] - n0["Ill"])
    cuttosize!(par.indices["ill"][new], n0["Ill"])
end

function cuttosize!(list, l)
    for _ in 1:(length(list)-l)
        pop!(list)
    end
end

function saveorpush!(list, index, value)
    try
        list[index] = value
    catch
        push!(list, value)
    end
end

function birth!(ps, par, indx, old, new)
    #choose two genetic configurations to mate
    rand!(par.rndm, par.indices["healthy"][old])
    dropzeros!(par.traits[old][par.rndm[1]])
    dropzeros!(par.traits[old][par.rndm[2]])
    #generate offsprings genetic configuration
    offspring!(indx, par, rand(par.MutationsPerBirth), old, new)
    #add to population size
    ps["PopSize"] += 1
    #add the individual to the current population state dictionary
    if DiploidModel2.ispropagable(par.traits[new][indx])
        #save index of new individual in appropriate list
        saveorpush!(par.indices["healthy"][new], ps["PopSize"] - ps["Ill"], indx)
    else
        ps["Ill"] += 1
        #save index of new individual in approriate list
        saveorpush!(par.indices["ill"][new], ps["Ill"], indx)
    end
    ps["ML"] += sum(par.traits[new][indx])
end

function offspring!(offspring_index, par, n_mut, old, new)
    #randomly recombine the parental genetic information
    par.traits[new][offspring_index] .=
        DiploidModel2.choose.(par.traits[old][par.rndm[1]], par.traits[old][par.rndm[2]])
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ = 1:n_mut
            mutate!(par, offspring_index, rand(par.MutationLocation), new)
        end
    end
    nothing
end

function mutate!(par, offspring_index, location, new)
    a = par.traits[new][offspring_index][location]
    #if there are still mutations to add
    if a < 2
        #add a mutation if a=0 and if a=1 add a mutation with prob. 1/2
        par.traits[new][offspring_index][location] = a * rand(Bool) + 1
    end
end

const_popsize(popsize,par,t) = par.K
random_fluctuations(popsize,par,t) = round(Integer,rand(par.fluctuations))

end # module WrightFisher

#---
