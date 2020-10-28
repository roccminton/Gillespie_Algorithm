using Random
using Distributions
"""
        run_gillespie(time::UnitRange{Int64},n₀::PopulationState,par::Parameter)

Run a exact stochastic simulation and t
"""
function run_gillespie(
    time::AbstractVector,
    n₀::PopulationState,
    par::Parameter,
    )

    #Setup
    current_time = time[1]
    current_population = n₀
    population_history = Vector{PopulationState}(undef,length(time))

    #run simulation
    for (index,step) in enumerate(time)
        while current_time ≤ step
            r = rates(current_population,par)
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
                    current_population = execute(i,current_population)
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
