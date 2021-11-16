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

export modelsetup,
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
function modelsetup(name::AbstractString,execute::Function)
    if name ∈ ["Yule","Linear","Logistic","Immigration",
        "YuleAcc","LinearAcc","LogisticAcc","ImmigrationAcc"]
        return (
            name = name,
            rates! = getfield(BirthDeath,Symbol(string(name,"Rates!"))),
            execute! = execute
            )
    else
        error("Unknown model name: $name")
    end
end

function modelsetup(name::AbstractString; rescaled::Bool=false)
    if rescaled
        modelsetup(string(name,"Acc"),execute_resc!)
    else
        modelsetup(name,execute_abs!)
    end
end

"Rate functions"
linear(n,bd) = n*bd
logisticdeath(n,d,c) = n*d+n*(n-1)*c
logisticbirth_rescaled(nK,d,c,K) = nK*(d+nK*c/K)
immigrationbirth(n,b,i) = n*b + i

"""
    LogisticRates(ps::pop_state,pr::par)

Return an array consisting of Logistic [birth,death] rates.
"""
function LogisticRates!(rates,n,pr)
    @inbounds rates[1] = linear(n,pr.birth)
    @inbounds rates[2] = logisticdeath(n,pr.death,pr.competition)
    nothing
end
"""
    LogisticAccRates(ps::pop_state,pr::par)

Return an array consisting of accelerated Logistic [birth,death] rates.
"""
function LogisticAccRates!(rates,n,pr)
    @fastmath n = n*pr.K
    @fastmath @inbounds rates[1] = linear(n,pr.birth)
    @fastmath @inbounds rates[2] = logisticdeath_rescaled(n,pr.death,pr.competition,pr.K)
    nothing
end

"""
    YuleRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
function YuleRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = linear(n,pr.birth)
    @fastmath @inbounds rates[2] = 0.0
    nothing
end
"""
    YuleAccRates(ps::pop_state,pr::par)

Return an array consisting of Yule [birth,death] rates.
"""
YuleAccRates!(rates,n,pr) = YuleRates!(rates,n*pr.K,pr)

"""
    LinearRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
function LinearRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = linear(n,pr.birth)
    @fastmath @inbounds rates[2] = linear(n,pr.death)
    nothing
end
"""
    LinearAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates.
"""
LinearAccRates!(rates,n,pr) = LinearRates!(rates,n*pr.K,pr)

"""
    ImmigrationRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
function ImmigrationRates!(rates,n,pr)
    @fastmath @inbounds rates[1] = immigrationbirth(n,pr.birth,pr.immigration)
    @fastmath @inbounds rates[2] = linear(n,pr.death)
    nothing
end

"""
    ImmigrationAccRates(ps::pop_state,pr::par)

Return an array consisting of Linear [birth,death] rates with immigration.
"""
ImmigrationAccRates!(rates,n,pr) = ImmigrationRates!(rates,n*pr.K,pr)

"""
    execute_diff(i,ps,diff,pr)

Executes the event that was choosen by the handed in integer 'i'.
Event 1 executes a birth event by adding 'diff' to the population state.
Event 2 executes a death event by decreasing the population state by 'diff'.
"""
function execute_diff!(i,ps,diff,pr)
    if i==1
        @fastmath n += diff
    elseif i==2
        @fastmath n -= diff
    else
        error("Index Error: No event #$i")
    end
    return n
end

execute_abs!(i,ps,pr) = execute_diff!(i,ps,1,pr)
execute_resc!(i,n,pr) = execute_diff!(i,ps,1/pr.K,pr)

end
