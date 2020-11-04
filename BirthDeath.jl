"""
    module BirthDeath

Defines all critical functions and types for a simple birth death process. Exported
types are
    Parameter -> Collection of parameters
    PopulationState -> Current population state
    ModelConfiguration -> essential model functions
Exported functions are
    Rates -> functions that return an array of rates for different events
    Execute ->  functions that execute the event (e.g. birth, death)
"""
module BirthDeath

using Parameters #for macro @with_kw

export Parameter, ModelConfiguration, PopulationState,
    execute_abs, execute_resc,
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

Immutable type of the current population state containing only the population size.
The population size can either be an integer, or a float if the size of individuals
is rescaled.
"""
struct PopulationState{T<:Real}
    "Population Size"
    size :: T
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
    if name in ["Yule","Linear","Logistic","Immigration",
        "YuleAcc","LinearAcc","LogisticAcc","ImmigrationAcc"]
        return ModelConfiguration(
            name,
            getfield(BirthDeath,Symbol(string(name,"Rates"))),
            execute
            )
    else
        error("Unknown model name: $name")
    end
end

function ModelConfiguration(name::AbstractString; rescaled::Bool=false)
    if rescaled
        ModelConfiguration(string(name,"Acc"),execute_resc)
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
    LogisticAccRates(ps::pop_state,pr::par)

Return an array consisting of accelerated Logistic [birth,death] rates.
"""
function LogisticAccRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size*pr.K
    return [n*pr.birth,n*(pr.death+n*pr.competition/pr.K)]
end
"""
    YuleRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    return [ps.size*pr.birth,0]
end
"""
    YuleAccRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleAccRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    return [ps.size*pr.birth*pr.K,0]
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
    LinearAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearAccRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size
    return [n*pr.birth*pr.K,n*pr.death*pr.K]
end
"""
    ImmigrationAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationAccRates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size
    return [n*pr.birth*K + pr.immigration,n*pr.death*K]
end

"""
    execute_abs(i::Int64,ps::PopulationState,pr::Parameter)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding one individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by one.
"""
function execute_abs(i::Int64,ps::PopulationState,pr::Parameter)::PopulationState
    if i==1
        return PopulationState{Int64}(ps.size + 1)
    elseif i==2
        return PopulationState{Int64}(ps.size - 1)
    else
        error("Index Error: No event #$i")
    end
end

"""
    execute_resc(i::Int64,ps::PopulationState,pr::Parameter)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding 1/K individuum to the population state size.
Event 2 executes a death event by decreasing the population state size by 1/K.
"""
function execute_resc(i::Int64,ps::PopulationState,pr::Parameter)::PopulationState
    if i==1
        return PopulationState{Float64}(ps.size + 1/pr.K)
    elseif i==2
        return PopulationState{Float64}(ps.size - 1/pr.K)
    else
        error("Index Error: No event #$i")
    end
end

end
