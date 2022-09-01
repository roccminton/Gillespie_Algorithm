"""

module DiploidModel2

Implementation of the same DiploidModel, but the genetic information from the
individuals in the population is saved only once in shared memory and then only
overwritten. Therefore 3 Arrays are implemented containing the indexes of the
individuals that are 1) alive and healthy 2) alive and ill and 3) not alive, hence
free space to use
"""

module DiploidModel2

using SparseArrays
using Random
using Distributions
using DataFrames

include("MainFunctions.jl")
import .Gillespie

function rungillespie(time,n₀,model_parameter)

    #setup empty population history
    l = length(time)
    population_history = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel2,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    #execute simulation
    Gillespie.run_gillespie!(
        time,n₀,
        setupparameter(model_parameter,n₀,l),
        execute!,
        rates!,
        initrates,
        population_history
        )

    return population_history

end

generatehealthypopulation(popsize) = initpopulation(popsize,0,0)

function initpopulation(popsize,ML,Ill)
    if ML >= 2*Ill
        return Dict(
            "PopSize" => popsize,
            "Ill" => Ill,
            "ML" => ML
            )
    else
        error("Mutation Load ($ML) must be at least twice the number of ill individual($Ill)")
    end
end

#---
#General Model Functions for parameter changes during runtime

function changeparameter(par,newpar)
    return merge(par,newpar)
end

function runseveralgillespie(totaltime,switch::Vector,parameter_changes,n₀,model_parameter)

    #setup empty population history
    l = length(totaltime)
    population_history = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel2,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    par = setupparameter(model_parameter,n₀,l)

    #save initial population
    for (x,nₓ) in n₀
        #add the population history for the given trait
        population_history[x][1] = nₓ
    end

    #switching times
    times = timesfromswitches(totaltime,switch)

    #execute simulation
    for (i,t) in enumerate(times)
        Gillespie.run_gillespie!(
            t,n₀,
            par,
            execute!,
            rates!,
            initrates,
            population_history,
            hstart = t[1]
            )
        #change parameter for next block, unless no block follows
        i < length(times) && (par = changeparameter(par,parameter_changes[i]))
    end
    return population_history
end

function runseveralgillespie(totaltime,switch::Number,par_changes,x0,parameters)
    return runseveralgillespie(totaltime,[switch],par_changes,x0,parameters)
end

#---
#Model Functions - Total Recombination

setupparameter(par,n0,historylength) = (
    par...,
    rndm = Vector{Int}(undef,2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = 1:par.Nloci,
    traits = inittraitsfromdict(par,n0),
    indices = Dict(
        "healthy" => collect(n0["Ill"]+1:n0["PopSize"]),
        "ill" => collect(1:n0["Ill"]),
        "free" => collect(n0["PopSize"]+1:round(Int,par.K + sqrt(par.K)))
    ),
    historylength = historylength
    )

function inittraitsfromdict(par,n0)
    locs = 1:par.Nloci
    traits = [spzeros(par.Nloci) for _ in 1:round(Int,par.K + sqrt(par.K))]
    for i in 1:n0["Ill"]
        l = rand(locs)
        traits[i][l] = 2
    end
    individuals = n0["Ill"]+1:n0["PopSize"]
    for i in n0["Ill"]+1:n0["ML"]-2*n0["Ill"]
        ind = rand(individuals)
        l = rand(locs)
        while traits[ind][l] ≠ 0
            ind = rand(individuals)
            l = rand(locs)
        end
        traits[ind][l] = 1
    end
    return traits
end

ispropagable(a::SparseVector) = !(2 ∈ a.nzval)

"""
    function rates!(rates,ps,par)

        > rates is a vector of the total birth and death rate as entries
        > ps is the population state as a sparse array
        > par is a named tuple with all relevant modle parameters
"""
function truerates!(rates, ps, par)
    #linear birth for all propagable individuals
    rates[1] = par.birth * (ps["PopSize"]-ps["Ill"])
    #uniform logistic death
    rates[2] = ps["PopSize"] * par.death + ps["PopSize"] * (ps["PopSize"] - 1) * par.competition
    nothing
end

function allbirthrates!(rates, ps, par)
    #linear birth for all individuals
    rates[1] = (par.birth * ps["PopSize"]) * (!isempty(par.indices["healthy"]))
    #uniform logistic death
    rates[2] = ps["PopSize"] * par.death + ps["PopSize"] * (ps["PopSize"] - 1) * par.competition
    nothing
end

function offspring!(offspring_index, par, n_mut)
    #randomly recombine the parental genetic information
    par.traits[offspring_index] .= choose.(
        par.traits[par.rndm[1]],
        par.traits[par.rndm[2]]
    )
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ in 1:n_mut
            mutate!(par,offspring_index,rand(par.MutationLocation))
        end
    end
    nothing
end

function mutate!(par,offspring_index,location)
    a = par.traits[offspring_index][location]
    #if there are still mutations to add
    if a < 2
        #add a mutation if a=0 and if a=1 add a mutation with prob. 1/2
        par.traits[offspring_index][location] = executemutation(a)
    end
end

function executemutation(existing_mutations)
    if iszero(existing_mutations)
        return 1
    else
        rand(Bool) ? (return 1) : (return 2)
    end
end

function choose(i, j)
    if i == j == 0
        return 0
    elseif iszero(i*j)
        return rand(Bool) ? 0 : 1
    else
        return (rand(Bool) ? 1 : (rand(Bool) ? 0 : 2))
    end
end

function birth!(ps, par)
    #choose two genetic configurations to mate
    rand!(par.rndm,par.indices["healthy"])
    dropzeros!(par.traits[par.rndm[1]]); dropzeros!(par.traits[par.rndm[2]])
    #select free index for offspring
    if isempty(par.indices["free"])
        offspring_index = length(par.traits) + 1
        push!(par.traits,spzeros(par.Nloci))
    else
        offspring_index = pop!(par.indices["free"])
    end
    #generate offsprings genetic configuration
    offspring!(offspring_index, par, rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    if ispropagable(par.traits[offspring_index])
        push!(par.indices["healthy"],offspring_index)
    else
        ps["Ill"] += 1
        push!(par.indices["ill"],offspring_index)
    end
    ps["PopSize"] += 1
    ps["ML"] += sum(par.traits[offspring_index])
end

function death!(ps, par)
    #choose fey
    if rand()<=ps["Ill"]/ps["PopSize"]
        fey_index = popat!(par.indices["ill"],rand(1:ps["Ill"]))
        ps["Ill"] -= 1
    else
        fey_index = popat!(par.indices["healthy"],rand(1:(ps["PopSize"]-ps["Ill"])))
    end
    #add fey to gravejard
    push!(par.indices["free"],fey_index)
    #update population state
    ps["PopSize"] -= 1
    ps["ML"] -= sum(par.traits[fey_index])
    nothing
end


function execute!(i,ps,par)
    if i == 1
        birth!(ps,par)
    elseif i == 2
        death!(ps,par)
    else
        error("Unknown event index: $i")
    end
end

isnotzero(n) = !iszero(n)

historytodataframe(t,history) = DataFrame(
    time=t,
    PopSize=history["PopSize"],
    Ill=history["Ill"],
    Mutation=history["ML"],
    )

function timesfromswitches(totaltime,switches)
    times = Vector{typeof(totaltime)}(undef,length(switches)+1)
    tstart = 1
    for (i,t) in enumerate(switches)
        tend = findfirst(x -> x≥t,totaltime)
        times[i] = totaltime[tstart:tend]
        tstart = tend+1
    end
    times[end] = totaltime[tstart:end]
    return times
end



end  # end of module DiploidModel2
