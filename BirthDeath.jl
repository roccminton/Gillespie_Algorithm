#=
    module BirthDeath

Defines all critical functions and types for a simple birth death process. Exported
types are
    Parameter -> Collection of parameters
    PopulationState -> Current population state
    PopulationStatistic -> Population State to save for statistics
    ModelConfiguration -> essential model functions
Exported functions are
    Rates -> functions that return an array of rates for different events
    Execute ->  functions that execute the event (e.g. birth, death)
=#

module BirthDeath

include("MainFunctions.jl")
import .Gillespie: ModelConfiguration

using Parameters #for macro @with_kw

export ModelConfiguration,
    execute_abs!, execute_resc!,
    LogisticRates, YuleRates, LinearRates, ImmigrationRates,
    LogisticAccRates, YuleAccRates, LinearAccRates, ImmigrationAccRates

#=
    ModelConfiguration(name::AbstractString; <keywords argument>)

Outer Constructor for the ModelConfiguration struc. Determines the rates and execute
Function from the name. Known model names are:
    Yule, Linear, Logistic, Immigration
Additionally a boolean keyword argument 'rescaled' can be handed in. If true the
`resc_execute` function is choosen, otherwise the `abs_executes`.
=#
function ModelConfiguration(name::AbstractString,execute::Function)
    if name ∈ ["Yule","Linear","Logistic","Immigration",
        "YuleAcc","LinearAcc","LogisticAcc","ImmigrationAcc"]
        return ModelConfiguration(
            name,
            getfield(BirthDeath,Symbol(string(name,"Rates!"))),
            execute
            )
    else
        error("Unknown model name: $name")
    end
end

function ModelConfiguration(name::AbstractString; rescaled::Bool=false)
    if rescaled
        ModelConfiguration(string(name,"Acc"),execute_resc!)
    else
        ModelConfiguration(name,execute_abs!)
    end
end

"""
    LogisticRates(ps::pop_state,pr::par)

Return an array consisting of Logistic [birth,death] rates.
"""
function LogisticRates!(rates,n,pr)
    @inbounds rates[1] = n*pr.birth
    @inbounds rates[2] = n*pr.death+n*(n-1)*pr.competition
    nothing
end
"""
    LogisticAccRates(ps::pop_state,pr::par)

Return an array consisting of accelerated Logistic [birth,death] rates.
"""
function LogisticAccRates!(rates,n,pr)
    @fastmath n = n*pr.K
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = n*(pr.death+n*pr.competition/pr.K)
    nothing
end

"""
    YuleRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = 0.0
    nothing
end
"""
    YuleAccRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleAccRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = n*pr.birth*pr.K
    @fastmath @inbounds rates[2] = 0.0
    nothing
end
"""
    LinearRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end
"""
    LinearAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearAccRates!(rates,n,pr)
    @fastmath n = n*pr.K
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end
"""
    ImmigrationRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = n*pr.birth + pr.immigration
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end

"""
    ImmigrationAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationAccRates!(rates,n,pr)
    @fastmath n = n*pr.K
    @fastmath @inbounds rates[1] = n*pr.birth + pr.immigration
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end

"""
    execute_abs(i::Int64,ps::PopulationState,pr::Parameter)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding one individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by one.
"""
function execute_abs!(i,ps,pr)
    if i==1
        @fastmath n += 1
    elseif i==2
        @fastmath n -= 1
    else
        error("Index Error: No event #$i")
    end
    return n
end

"""
    execute_resc(i::Int64,ps::PopulationState,pr::Parameter)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding 1/K individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by 1/K.
"""
function execute_resc!(i,n,pr)
    if i==1
        @fastmath n += 1/pr.K
    elseif i==2
        @fastmath n -= 1/pr.K
    else
        error("Index Error: No event #$i")
    end
    return n
end

end
