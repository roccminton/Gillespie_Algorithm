"""
Implementation of a simple Wright Fisher model with constant
population size
"""

module WrightFisher

using Random
using Distributions
using ProgressMeter
using SparseArrays
using Random
using SpecialFunctions

#execute simulation
function run_wrightfisher(time, n₀, model_parameter)

    #setup empty population history
    l = length(time)
    population_history = setup_pop_hist(model_parameter,n₀,l)
    parameter = setupparameter(model_parameter, n₀, l)

    #execute simulation
    mainiteration!(
        n₀,
        time,
        parameter,
        population_history,
    )

    return population_history
end



#placeholder to be overwritten in specification
function setupparameter end
function birth! end
function setup_pop_hist end
function saveonestep! end

function mainiteration!(n0, time, par, pop_hist)
    #run simulation
    @showprogress for t in time
        #save one step evolution
        saveonestep!(pop_hist, t, n0, par)
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
    for (k,v) in n0
        cuttosize!(par.indices[k][new], v)
    end
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

end # module WrightFisher

#---
