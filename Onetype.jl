#=
    module OneType

Defines all critical functions and types for a simple birth death process.

Rates are represended as Vector with two entries only, where the first
entry is the overall birht rate and the second the death rate of the system.

=#

module OneType

include("MainFunctions.jl")
import .Gillespie

export modelsetup, rungillespie,
    execute_abs, execute_resc,
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
            rates! = getfield(OneType,Symbol(string(name,"Rates!"))),
            execute! = execute
            )
    else
        error("Unknown model name: $name")
    end
end

function modelsetup(name::AbstractString; rescaled::Bool=false)
    if rescaled
        modelsetup(string(name,"Acc"),execute_resc)
    else
        modelsetup(name,execute_abs)
    end
end

"Rate functions"
linear(n,bd) = n*bd
logisticdeath(n,d,c) = n*d+n*(n-1)*c
logisticdeath_rescaled(nK,d,c,K) = nK*(d+nK*c/K)
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
function execute_diff(i,ps,diff,pr)
    if i==1
        @fastmath ps += diff
    elseif i==2
        @fastmath ps -= diff
    else
        error("Index Error: No event #$i")
    end
    return ps
end

execute_abs(i,ps,pr) = execute_diff(i,ps,1,pr)
execute_resc(i,ps,pr) = execute_diff(i,ps,1/pr.K,pr)

setuprates(birthrate) = Vector{typeof(birthrate)}(undef,2)
setuphistory(time,n₀) = zeros(typeof(n₀),length(time))

function rungillespie(t,n_0,par,conf)
    #setup empty history
    population_history = setuphistory(t,n_0)
    #execute simulation
    Gillespie.run_gillespie!(
        t,n_0,par,
        conf.execute!,conf.rates!,
        setuprates(par.birth),
        population_history
        )

    return population_history
end

end #end of module
