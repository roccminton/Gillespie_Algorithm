module Gillespie

using Random
using Distributions
using ProgressMeter

"""
        run_gillespie(time::AbstractVector,n₀,par::Parameter)

Run a exact stochastic simulation over the Interval time with an initial Population
`n₀`, Parameter Set `par` and Model Configuration `conf`. The output of the function
is a Vector{PopulationState} of length `length(time)` with the history of the
states of the population during the simulation.
"""
function run_gillespie!(time,n₀,par,execute!,rates!,initrates,population_history)

    mainiteration!(
        population_history,
        initrates,
        n₀,
        convert(Float64,time[1]),
        time,
        par,
        execute!,
        rates!
    )
    
end

function mainiteration!(pop_hist,rates,n0::Real,ct,time,par,ex!,r!)
    #run simulation
    @showprogress for (index,step) in enumerate(time)
        #save one step evolution
        pop_hist[index] = n0
        #execute one step of the simulation
        ct, n0 = onestep!(n0,rates,ct,step,par,ex!,r!)
        #check if step was completed or evolution stopped inbetween
        @fastmath step-ct > 0.0 && break
    end
end

function mainiteration!(pop_hist,rates,n0::Dict,ct,time,par,ex!,r!)
    #run simulation
    @showprogress for (index,step) in enumerate(time)
        #save one step evolution
        for (x,vₓ) in n0
            !haskey(pop_hist,x) && (pop_hist[x] = zeros(eltype(valtype(pop_hist)),length(time)))
            pop_hist[x][index] = vₓ[1]
        end
        #execute one step of the simulation
        ct = onestep!(n0,rates,ct,step,par,ex!,r!)
        #check if step was completed or evolution stopped inbetween
        @fastmath step-ct > 0.0 && break
    end
end

function mainiteration!(pop_hist,rates,n0,ct,time,par,ex!,r!)
    #run simulation
    @showprogress for (index,step) in enumerate(time)
        #save one step evolution
        view(pop_hist,index,:) .= n0
        #execute one step of the simulation
        ct = onestep!(n0,rates,ct,step,par,ex!,r!)
        #check if step was completed or evolution stopped inbetween
        @fastmath step-ct > 0.0 && break
    end
end

"""
Executes one step of the evolution by modifying `x_0` and `rates`.
"""
function onestep!(x_0,rates::Matrix,t_0,t_end,par,ex!,r!)
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        ex!(i,x_0,par)
    end
    return t_0
end

function onestep!(x_0,rates::Dict,t_0,t_end,par,ex!,r!)
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, trait, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        ex!(i,trait,x_0,rates,par)
    end
    return t_0
end

function onestep!(x_0::Real,rates::Vector,t_0,t_end,par,ex!,r!)
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        x_0 = ex!(i,x_0,par)
    end
    return t_0, x_0
end

"""
    nexteventandtime(rates::Vector{Float64})

Samples a exponential distributed random variable to determine the time for the next
event and calls `choose_event`. The return value is a tuple consiting of the
envent index returned by `choose_event` and the time to the next event.
"""
function nexteventandtime(rates)
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

function nexteventandtime(rates::Dict)
    #calculate total event rate
    total_rate = sumsumdict(rates)
    #Population in absorbing state if sum of rates is zero
    iszero(total_rate) && return 0, zero(keytype(rates)), 0.0
    #sample next event time
    dt = rand(Exponential(1 / (total_rate)))
    #choose event
    i, trait = chooseevent(rates,total_rate)
    return i, trait, dt
end

"""
    chooseevent(rates::Vector{Float64},total_rate::Float64;<keyword arguments>)

From the vector of total rates choose at random one of the indices of the vector
according to their rates.
The value 0 is returned if the total rates are positive, but too smale to let
the evolution continue. The maximum number of tries is set by `max_try=1000`.
"""
function chooseevent(rates,total_rate)
    #make it a uniform random variable in (0,total_rate)
    rndm = rand(Uniform(0.0,total_rate))
    #choose the rate at random
    @inbounds for (index, rate) in enumerate(rates)
        rndm -= rate
        rndm ≤ 0.0 && return index
    end
end

function chooseevent(rates::Dict,total_rate)
    #make it a uniform random variable in (0,total_rate)
    rndm = rand(Uniform(0.0,total_rate))
    #choose the rate at random
    for (trait, rate) in rates
        for (index,r) in enumerate(rate)
            rndm -= r
            rndm ≤ 0.0 && return index, trait
        end
    end
end

"""
For a dictionary with vectors of values calculates the sum of all
values of all vectors combined.
"""
function sumsumdict(D)
    result = zero(eltype(valtype(D)))
    for v in values(D)
        result += sum(v)
    end
    return result
end

end #end of module Gillespie
