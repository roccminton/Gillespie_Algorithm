"""
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
"""
module BirthDeath

using Parameters #for macro @with_kw

export Parameter, ModelConfiguration, PopulationState, PopulationStatistic,
    generate_statistic,
    execute_abs!, execute_resc!,
    LogisticRates, YuleRates, LinearRates, ImmigrationRates,
    LogisticAccRates, YuleAccRates, LinearAccRates, ImmigrationAccRates

"""
    struct Parameter

Compilation of immutable parameters for the whole runtime of the simulation.
"""
@with_kw struct Parameter{T<:AbstractFloat}
    "birth rate"
    birth :: T = 0.0
    "death rate"
    death :: T = 0.0
    "competition preassure"
    competition :: T = 0.0
    "carrying capacity"
    K :: T = 1.0
    "immigration rate"
    migration :: T = 0.0
end

"""
    struct PopulationState

Mutable type of the current population state containing only the population size.
The population size can either be an integer, or a float if the size of individuals
is rescaled.
"""
mutable struct PopulationState{T<:Real}
    "Population Size"
    size :: T
end

struct PopulationStatistic{T<:Real}
    "Population Size"
    size :: T
end

function PopulationStatistic(ps::PopulationState{T}) where T
    return PopulationStatistic{T}(ps.size)
end

"""
    struc ModelConfiguration

Declares configurations of the model that unlike the model parameters don't have
to be passed to the rates function, but setup the model initially.
"""
struct ModelConfiguration
    name :: AbstractString
    rates! :: Function
    execute! :: Function
end

"""
    ModelConfiguration(name::AbstractString; <keywords argument>)

Outer Constructor for the ModelConfiguration struc. Determines the rates and execute
Function from the name. Known model names are:
    Yule, Linear, Logistic, Immigration
Additionally a boolean keyword argument 'rescaled' can be handed in. If true the
`resc_execute` function is choosen, otherwise the `abs_executes`.
"""
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
function LogisticRates!(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    n = ps.size
    @inbounds rates[1] = n*pr.birth
    @inbounds rates[2] = n*pr.death+n*(n-1)*pr.competition
    nothing
end
"""
    LogisticAccRates(ps::pop_state,pr::par)

Return an array consisting of accelerated Logistic [birth,death] rates.
"""
function LogisticAccRates!(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath n = ps.size*pr.K
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = n*(pr.death+n*pr.competition/pr.K)
    nothing
end

"""
    YuleRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleRates(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath @inbounds rates[1] = ps.size*pr.birth
    @fastmath @inbounds rates[2] = 0.0
    nothing
end
"""
    YuleAccRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleAccRates(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath @inbounds rates[1] = ps.size*pr.birth*pr.K
    @fastmath @inbounds rates[2] = 0.0
    nothing
end
"""
    LinearRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearRates(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath n = ps.size
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end
"""
    LinearAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearAccRates(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath n = ps.size*pr.K
    @fastmath @inbounds rates[1] = n*pr.birth
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end
"""
    ImmigrationRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationRates(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath n = ps.size
    @fastmath @inbounds rates[1] = n*pr.birth + pr.immigration
    @fastmath @inbounds rates[2] = n*pr.death
    nothing
end

"""
    ImmigrationAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationAccRates(rates::Vector{Float64},ps::PopulationState,pr::Parameter)
    @fastmath n = ps.size*pr.k
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
function execute_abs!(i::Int64,ps::PopulationState,pr::Parameter)
    if i==1
        @fastmath ps.size += 1
    elseif i==2
        @fastmath ps.size -= 1
    else
        error("Index Error: No event #$i")
    end
    nothing
end

"""
    execute_resc(i::Int64,ps::PopulationState,pr::Parameter)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding 1/K individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by 1/K.
"""
function execute_resc!(i::Int64,ps::PopulationState,pr::Parameter)
    if i==1
        @fastmath ps.size += 1/pr.K
    elseif i==2
        @fastmath ps.size -= 1/pr.K
    else
        error("Index Error: No event #$i")
    end
    nothing
end

end
