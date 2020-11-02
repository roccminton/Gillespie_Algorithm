"""
    module BirthDeath

Defines all critical functions and types for a simple birth death process. Exported
types are
    par -> Collection of parameters
    pop_state -> Current population state
Exported functions are
    rates(pop_state,par) -> returns an array of rates for different events
    events -> array of functions that execute the event (e.g. birth, death)
"""
module BirthDeath

using Parameters

export Parameter, ModelConfiguration, PopulationState,
    execute,
    LogisticRates, YuleRates, LinearRates, ImmigrationRates

"""
    struct Parameter

Compilation of immutable parameters for the whole runtime of the simulation.
"""
@with_kw struct Parameter
    "birth rate"
    birth :: Float64 = 0.0
    "death rate"
    death :: Float64 = 0.0
    "competition preassure"
    competition :: Float64 = 0.0
    "carrying capacity"
    K :: Int64 = 1
    "immigration rate"
    migration :: Float64 = 0.0
end

"""
    struct PopulationState

Immutable type of the current population state containing only the population size.
The population size can either be an integer, or a float if the size of individuals
is rescaled.
"""
struct PopulationState
    "Population Size"
    size :: Int64
end

"""
    struc ModelConfiguration

Declares configurations of the model that unlike the model parameters don't have
to be passed to the rates function, but setup the model initially.
"""
struct ModelConfiguration
    name :: AbstractString
    rates :: Function
    execute :: Function
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
    if name == "Yule"
        return ModelConfiguration(name,YuleRates,execute)
    elseif name == "Linear"
        return ModelConfiguration(name,LinearRates,execute)
    elseif name == "Logistic"
        return ModelConfiguration(name,LogisticRates,execute)
    elseif name == "Immigration"
        return ModelConfiguration(name,ImmigrationRates,execute)
    else
        error("Unknown model name: $name")
    end
end

function ModelConfiguration(name::AbstractString; rescaled::Bool=false)
    if rescaled
        ModelConfiguration(name,execute_resc)
    else
        ModelConfiguration(name,execute_abs)
    end
end

"""
    LogisticRates(ps::pop_state,pr::par)

Return an array consisting of Logistic [birth,death] rates.
"""
function LogisticRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size
    return [n*pr.birth,n*pr.death+n*(n-1)*pr.competition]
end
"""
    YuleRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    return [ps.size*pr.birth,0]
end
"""
    LinearRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size
    return [n*pr.birth,n*pr.death]
end
"""
    ImmigrationRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size
    return [n*pr.birth + pr.immigration,n*pr.death]
end

"""
    execute_abs(i::Int64,ps::pop_state)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding one individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by one.
"""
function execute_abs(i::Int64,ps::PopulationState)::PopulationState
    if i==1
        return PopulationState(ps.size + 1)
    elseif i==2
        return PopulationState(ps.size - 1)
    else
        error("Index Error: No event #$i")
    end
end

"""
    execute_resc(i::Int64,ps::pop_state)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding one individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by one.
"""
function execute_resc(i::Int64,ps::PopulationState,pr::Parameter)::PopulationState
    if i==1
        return PopulationState(ps.size + 1/pr.K)
    elseif i==2
        return PopulationState(ps.size - 1/pr.K)
    else
        error("Index Error: No event #$i")
    end
end

end
