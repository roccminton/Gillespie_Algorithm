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

export Parameter, PopulationState, execute, rates

"""
    struct Parameter

Compilation of immutable parameters for the whole runtime of the simulation.
"""
struct Parameter
    "birth rate"
    b :: Float64
    "death rate"
    d :: Float64
    "competition preassure"
    c :: Float64
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
    rates(ps::pop_state,pr::par)

Return an array consisting of [birth,death] rates.
"""
function rates(ps::PopulationState,pr::Parameter) :: Vector{Float64}
    n = ps.size
    return [n*pr.b,n*pr.d+n*(n-1)*pr.c]
end
"""
    birth(ps::pop_state)

Return a new pop_state instance with an increased population size by one
individuum.
"""
function birth(ps::PopulationState)::PopulationState
    return PopulationState(ps.size + 1)
end

"""
    death(ps::pop_state)

Return a new pop_state instance with an decreased population size by one
individuum.
"""
function death(ps::PopulationState)::PopulationState
    return PopulationState(ps.size - 1)
end
"""
    execute(i::Int64,ps::pop_state)

Collects all defined events for this model in the same order as the rates vector.
"""
function execute(i::Int64,ps::PopulationState)::PopulationState
    if i==1
        birth(ps)
    elseif i==2
        death(ps)
    else
        error("Index Error: No event #$i")
    end
end


end
