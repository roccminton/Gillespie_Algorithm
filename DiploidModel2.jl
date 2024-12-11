"""

module DiploidModel2

Implementation of the same DiploidModel, but the genetic information from the
individuals in the population is saved only once in shared memory and then only
overwritten. Therefore 3 Arrays are implemented containing the indexes of the
individuals that are 1) alive and healthy 2) alive and ill and 3) not alive, hence
free space to use

birth = b,
death = d,
competition = (b-d)/K,
μ = dni,
Nloci = N,
K = K,
rates = "allbirthrates!"
"""

module DiploidModel2

using SparseArrays
using Random
using Distributions
using DataFrames
using StatsBase

include("MainFunctions.jl")
import .Gillespie


"""
    Executes one run of the Gillespie algorithm with initial population size n0 over
    the time intervall time. The named tuple model_parameter has to have at least the following
    keys: birth, death, competition, μ, Nloci, K, recombination
"""
function rungillespie(time,n₀,model_parameter)

    l = length(time)

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel2,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    #add additional parameter for backend
    parameter = setupparameter(model_parameter,n₀,l)

    #add additional info to the parameter tuple according to the history type
    parameter = addstatsparameter(parameter,n₀,l)

    #setup empty population history
    population_history = setup_pop_hist(parameter,n₀,l)

    #choose statistic function
    statistic! = choosestatsfunction(population_history)

    #execute simulation
    Gillespie.run_gillespie!(
        time,n₀,
        parameter,
        execute!,
        rates!,
        initrates,
        population_history,
        statistic! = statistic!
        )

    #close any open files
    haskey(parameter,:file) && close(parameter.file)

    return population_history
end

#---

function setupparameter end
function setup_pop_hist end

#use the default stats function from Gillespie
choosestatsfunction(population_history) = Gillespie.saveonestep!
#nothing to add
addstatsparameter(ph,par,n0,l) = par

#--- Rates Function
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

#--- Execute Functions

function birth! end

function updatestats_death! end
function updatestats_birth! end

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
    updateps_death!(ps,par,fey_index)
end

function updateps_birth!(ps,par,offspring_index)
    if ispropagable(par.traits[offspring_index],par.Nloci)
        push!(par.indices["healthy"],offspring_index)
    else
        ps["Ill"] += 1
        push!(par.indices["ill"],offspring_index)
    end
    ps["PopSize"] += 1
    ps["ML"] += mutationload(par.traits[offspring_index])
    #update further statistics
    updatestats_birth!(ps,par,offspring_index)
end

function updateps_death!(ps,par,fey_index)
    ps["PopSize"] -= 1
    ps["ML"] -= mutationload(par.traits[fey_index])
    #update further statistics
    updatestats_death!(ps,par,fey_index)
end

ispropagable(a::SparseVector{Int64,Int64}) = !(2 ∈ a.nzval)
ispropagable(a::SparseVector{Int64,Int64},Nloci) = ispropagable(a)
function ispropagable(a::SparseVector{Bool,Int64},Nloci)
    for p in 1:Nloci
        a[p] && a[Nloci+p] && return false
    end
    return true
end
ispropagable(a::Vector,Nloci) = ispropagable(a)
function ispropagable(a::Vector)
    for (i,p) in enumerate(a[1])
        isone(p) && isone(a[2][i]) && return false
    end
    return true
end

mutationload(a::SparseVector) = sum(a)
mutationload(a::Vector) = sum(sum(svec) for svec in a)

#---

function execute!(i,ps,par)
    if i == 1
        birth!(ps,par)
    elseif i == 2
        death!(ps,par)
    else
        error("Unknown event index: $i")
    end
end

end  # end of module DiploidModel2
