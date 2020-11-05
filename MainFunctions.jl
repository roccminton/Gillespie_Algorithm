using Random
using Distributions
using ProgressMeter
"""
        run_gillespie(time::AbstractVector,n₀::PopulationState,par::Parameter)

Run a exact stochastic simulation over the Interval time with an initial Population
`n₀`, Parameter Set `par` and Model Configuration `conf`. The output of the function
is a Vector{PopulationState} of length `length(time)` with the history of the
states of the population during the simulation.
"""
function run_gillespie(
    time::UnitRange{Int64},
    n₀::PopulationState,
    par::Parameter,
    conf::ModelConfiguration
    )

    #Setup
    current_time = convert(Float64,time[1])
    current_population = n₀
    population_history = Vector{PopulationState}(undef,length(time))

    #run simulation
    return run_body(
        time,
        current_time,
        current_population,
        population_history,
        par,
        conf
        )
end

function run_body(time::UnitRange{Int64},
    t_0::Float64,
    n_0::PopulationState,
    history::Vector{PopulationState},
    par::Parameter,
    conf::ModelConfiguration
    )
    @showprogress for(index,step) in enumerate(time)
        t_0, n_0 = one_step(t_0,step,n_0,par,conf)
        history[index] = n_0
        step-t_0 > 0 && break
    end
    return history
end

function one_step(
    t_0 :: Float64,
    t_end :: Int64,
    x_0 :: PopulationState,
    par :: Parameter,
    conf :: ModelConfiguration
    )
    while t_0 ≤ t_end
        r = conf.rates(x_0,par)
        total_rate = sum(r)
        #Population in absorbing state if sum of rates is zero
        iszero(total_rate) && break
        #sample next event time
        dt = rand(Exponential(1 / (total_rate)))
        #update time
        t_0 += dt
        #choose event
        i = choose_event(r,total_rate)
        #execute event
        x_0 = conf.execute(i,x_0,par)
    end
    return t_0, x_0
end

"""
    choose_event(rates::Vector{Float64},total_rate::Float64;<keyword arguments>)

From the vector of total rates choose at random one of the indices of the vector
according to their rates.
The value 0 is returned if the total rates are positive, but too smale to let
the evolution continue. The maximum number of tries is set by `max_try=1000`.
"""
function choose_event(rates::Vector{Float64},total_rate::Float64;max_try::Int64=1000)::Int64
    #make it a uniform random variable in (0,total_rate)
    rndm = rand(Uniform(0.0,total_rate))
    #choose the rate at random
    for (index, rate) in enumerate(rates)
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
    n₀::PopulationState,
    pars::Vector{Parameter},
    confs::Vector{ModelConfiguration}
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
    n₀::PopulationState,
    pars::Vector{Parameter},
    conf::ModelConfiguration
    ) = run_gillespie(times,n₀,pars,fill(conf,length(times)))
