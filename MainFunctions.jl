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
    current_time = time[1]
    current_population = n₀
    population_history = Vector{typeof(n₀)}(undef,length(time))

    #run simulation
    @showprogress for (index,step) in enumerate(time)
        while current_time ≤ step
            r = conf.rates(current_population,par)
            total_rate = sum(r)
            #Population in absorbing state if sum of rates is zero
            if total_rate > 0
                #sample next event time
                dt = rand(Exponential(1 / (total_rate)))
                #update time
                current_time += dt
                #choose event
                i = choose_event(r,total_rate)
                if i > 0
                    current_population = conf.execute(i,current_population,par)
                else
                    break
                end
            else
                break
            end
        end
        #safe population status
        population_history[index] = current_population
    end
    return population_history
end

"""
    choose_event(rates::Vector{Float64},total_rate::Float64;<keyword arguments>)

From the vector of total rates choose at random one of the indices of the vector
according to their rates.
The value 0 is returned if the total rates are positive, but too smale to let
the evolution continue. The maximum number of tries is set by `max_try=1000`.
"""
function choose_event(rates::Vector{Float64},total_rate::Float64;max_try::Int64=1000)::Int64
    #setup a uniform random variable in (0,1)
    rndm = rand()
    #make it a uniform random variable in (0,total_rate)
    if total_rate ≥ 1
        #or by multiplying with the total_rate if its greater than 1
        rndm = rndm*total_rate
    else
        #or by rejecting samples if they are to big
        #set the maximum numbers of trys to 1000
        #otherwise stop the evolution
        broke = false
        for _ in 1:max_try
            if rndm < total_rate
                broke = true
                break
            else
                rndm = rand()
            end
        end
        if broke == false
            return 0
        end
    end
    for (index, rate) in enumerate(rates)
        rndm -= rate
        if rndm ≤ 0
            return index
        end
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
