module Gillespie

using Random
using Distributions
using ProgressMeter

export run_gillespie

"""
        run_gillespie(time::AbstractVector,n₀,par::Parameter)

Run a exact stochastic simulation over the Interval time with an initial Population
`n₀`, Parameter Set `par` and Model Configuration `conf`. The output of the function
is a Vector{PopulationState} of length `length(time)` with the history of the
states of the population during the simulation.
"""
function run_gillespie(time,n₀,par,conf)
    #Setup
    ct = convert(Float64,time[1]) #current time
    population_history = Vector{typeof(n₀)}(undef,length(time))
    rates = Vector{Float64}(undef,2)
    ex! = conf.execute!
    r! = conf.rates!

    if typeof(n₀).mutable
        #run simulation
        @showprogress for (index,step) in enumerate(time)
            #save one step evolution
            @inbounds population_history[index] = copy(n₀)
            #execute one step of the simulation
            ct = onestep_mutable!(n₀,rates,ct,step,par,ex!,r!)
            #check if step was completed or evolution stopped inbetween
            @fastmath step-ct > 0.0 && break
        end
    else
        #run simulation
        @showprogress for (index,step) in enumerate(time)
            #save one step evolution
            @inbounds population_history[index] = n₀
            #execute one step of the simulation
            ct, n₀ = onestep_immutable!(n₀,rates,ct,step,par,ex!,r!)
            #check if step was completed or evolution stopped inbetween
            @fastmath step-ct > 0.0 && break
        end
    end
    return population_history
end

"""
    onestep!(
        t_0 :: Float64,
        t_end :: Int64,
        x_0 :: PopulationState,
        rates :: Vector{Float64},
        par :: Parameter,
        execute! :: Function,
        rates! :: Function
        )

Executes one step of the evolution by modifying `x_0` and `rates`.
"""
function onestep_mutable!(
    x_0 ,
    rates :: Vector{Float64},
    t_0 :: Float64,
    t_end :: Int64,
    par,
    execute! :: Function,
    rates! :: Function
    )

    while t_0 ≤ t_end
        rates!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        execute!(i,x_0,par)
    end
    return t_0
end

function onestep_immutable!(
    x_0 ,
    rates,
    t_0,
    t_end,
    par,
    execute!,
    rates!
    )

    while t_0 ≤ t_end
        rates!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        x_0 = execute!(i,x_0,par)
    end
    return t_0, x_0
end

"""
    nexteventandtime(rates::Vector{Float64})

Samples a exponential distributed random variable to determine the time for the next
event and calls `choose_event`. The return value is a tuple consiting of the
envent index returned by `choose_event` and the time to the next event.
"""
function nexteventandtime(rates::Vector{Float64})
    #calculate total event rate
    @fastmath total_rate = sum(rates)
    #Population in absorbing state if sum of rates is zero
    iszero(total_rate) && return 0, 0.0
    #sample next event time
    dt = rand(Exponential(1 / (total_rate)))
    #choose event
    i = chooseevent(rates,total_rate)
    return i, dt
end

"""
    chooseevent(rates::Vector{Float64},total_rate::Float64;<keyword arguments>)

From the vector of total rates choose at random one of the indices of the vector
according to their rates.
The value 0 is returned if the total rates are positive, but too smale to let
the evolution continue. The maximum number of tries is set by `max_try=1000`.
"""
function chooseevent(rates::Vector{Float64},total_rate::Float64)::Int64
    #make it a uniform random variable in (0,total_rate)
    rndm = rand(Uniform(0.0,total_rate))
    #choose the rate at random
    @inbounds for (index, rate) in enumerate(rates)
        rndm -= rate
        rndm ≤ 0.0 && return index
    end
end

"""
    run_gillespie(time::Vector{AbstractVector},n₀::PopulationState,par::Vector{Parameter},conf::Vector{ModelConfiguration})

Execute several 'runn_gillespie' calls in a row. Starting with the handed in initial
Population `n₀` and the first index in the `par` and `conf` Vector. After the first
`time` has expired a new simulation starts with the last `PopulationState` as initial
Population for the next run, and the second indices of the Vector etc. The return
is one long population history of all the runs concatenated together.
"""
function run_gillespie(
    times::Vector{UnitRange{Int64}},
    n₀,
    pars,
    confs
    )

    #Check if all Vector have the same length
    if length(times)==length(pars)==length(confs)
        #Setup
        population_history = Vector{PopulationState}(undef,0)
        sub_history = Vector{PopulationState}([n_0])
        #run simulations in row
        for (i, time) in enumerate(times)
            sub_history = run_gillespie(
                                time,
                                sub_history[end],
                                pars[i],
                                confs[i]
                                )
            #add subhistorys
            population_history = vcat(population_history,sub_history)
        end
        return population_history
    else
        error(
            """Handed in Vectors in run_gillespie have to be of equal length.
            time: $(length(times)), par: $(length(pars)), conf: $(length(confs))."""
            )
    end
end

"""
    run_gillespie(times::Vector{AbstractVector},n₀::PopulationState,pars::Vector{Parameter},conf::ModelConfiguration)

Assume that if the model Configuration doesn't come in an Vector the same
Configuration should be used for all time intervals.
"""
run_gillespie(
    times::Vector{UnitRange{Int64}},
    n₀,
    pars,
    conf
    ) = run_gillespie(times,n₀,pars,fill(conf,length(times)))

end
