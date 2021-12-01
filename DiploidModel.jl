"""

module DiploidModel

Implementation of a diploid model, where individuals are characterized
with a number *N* of mutable, diploid genes resp. gene sections.
Mutations are one way only, recessive and leathal, meaning that
any individual carrying two mutations at one position is excluded
from the reproduction process. Besides that every individual has
the same birth and death rate.

Internally we represent the type of an individual as a struct holding the genetic
information at the field *genes* as a sparse Vector of length N

"""

module DiploidModel

using SparseArrays
using Distributions

include("MainFunctions.jl")
import .Gillespie

function rungillespie(time,n₀,model_parameter)

    #setup empty population history
    l = length(time)
    population_history = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #execute simulation
    Gillespie.run_gillespie!(
        time,n₀,
        setupparameter(model_parameter,n₀),
        execute!,
        rates!,
        initrates,
        population_history
        )

    return population_history

end

generatehealthypopulation(popsize,Nloci) = Dict(RandomIndividual(Nloci) => popsize)


setupparameter(par,n0) = (
    par...,
    rndm = rand(2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = DiscreteUniform(1,par.Nloci),
    parent1 = spzeros(Int8,par.Nloci),
    parent2 = spzeros(Int8,par.Nloci),
    offspring = spzeros(Int8,par.Nloci),
    popsize = [sum(values(n0)),npropagable(n0)]
    )


struct RandomIndividual{T <: Integer}
    genes :: SparseVector{Int8,T}
    ispropagable :: Bool
end

RandomIndividual(a::SparseVector) = RandomIndividual(
    sparsevec(a.nzind,Int8.(a.nzval),length(a)),
    ispropagable(a)
    )

RandomIndividual(N::Integer) = RandomIndividual(
    spzeros(Int8,N),
    true
)

ispropagable(a::SparseVector) = !(2 ∈ a.nzval)

function isindictandadd!(d,x,diff)
    isindict = false
    for y in keys(d)
        if y.genes == x
            d[y] += diff
            isindict = true
            break
        end
    end
    !isindict && (d[RandomIndividual(x)] = diff)
end

"""
    Sum of all healty individual in the population state ps
"""
npropagable(ps) = sum(x.ispropagable && nₓ for (x,nₓ) in ps)

"""
    function rates!(rates,ps,par)

        > rates is a vector of the total birth and death rate as entries
        > ps is the population state as a sparse array
        > par is a named tuple with all relevant modle parameters
"""
function rates!(rates, ps, par)
    #calculate total population size
    par.popsize[1] = sum(values(ps))
    #calculate size of propagable subpopulation
    par.popsize[2] = npropagable(ps)
    #linear birth for all propagable individuals
    rates[1] = par.birth * par.popsize[2]
    #uniform logistic death
    rates[2] = par.popsize[1] * par.death + par.popsize[1] * (par.popsize[1] - 1) * par.competition
    nothing
end

function choosecouple!(rndm, ps, par)
    #make it a uniform random variable in (0,total_rate)
    U = Uniform(0.0, par.popsize[2])
    rndm[1] = rand(U)
    rndm[2] = rand(U)
    #initiate
    searchparent1 = true
    searchparent2 = true
    #choose the rate at random
    for (x, nₓ) in ps
        x.ispropagable && (rndm .-= nₓ)
        searchparent1 && rndm[1] ≤ 0.0 && (par.parent1 .= x.genes; searchparent1 = false)
        searchparent2 && rndm[2] ≤ 0.0 && (par.parent2 .= x.genes; searchparent2 = false)
        !searchparent1 && !searchparent2 && break
    end
    nothing
end

function choosefey(ps,pop_size)
    #make it a uniform random variable in (0,total_rate)
    rndm = rand(Uniform(0.0,pop_size))
    #choose the rate at random
    for (x, nₓ) in ps
        rndm -= nₓ
        rndm ≤ 0.0 && return x
    end
end

function offspring!(par, n_mut)
    #randomly recombine the parental genetic information
    par.offspring .= choose.(par.parent1, par.parent2)
    #add n_mut mutations to random positions mutation
    for _ = 1:n_mut
        mutate!(par, rand(par.MutationLocation))
    end
    nothing
end

function mutate!(par, location)
    #add a mutation if no mutation is present
    if par.offspring[location] == 0
        par.offspring[location] = 1
        #if there is already a mutation at the location choose at random if the same
        #copy mutates again, hence nothing happens or the other site mutates
    elseif par.offspring[location] == 1
        par.offspring[location] = rand(Bool) ? 1 : 2
    end
end

function choose(i, j)
    if i == j == 0
        return 0
    elseif (i == 0 && j == 1) || (i == 1 && j == 0)
        return rand(Bool) ? 0 : 1
    elseif i == j == 1
        return (rand(Bool) ? 1 : (rand(Bool) ? 0 : 2))
    else
        error("Cannot couple genetic information $i and $j")
    end
end

function birth!(ps, par)
    #choose two genetic configurations to mate
    choosecouple!(par.rndm, ps, par)
    dropzeros!(par.parent1); dropzeros!(par.parent2)
    #generate offsprings genetic configuration
    offspring!(par,rand(par.MutationsPerBirth))
    isindictandadd!(ps,par.offspring,one(valtype(ps)))
end

function death!(ps, par)
    ps[choosefey(ps,par.popsize[1])] -= one(valtype(ps))
end

function execute!(i,ps,par)
    if i == 1
        birth!(ps,par)
    elseif i == 2
        death!(ps,par)
    else
        error("Unknown event index: $i")
    end
end

"""
    populationsize

    Get the population size process from the history dictionary
"""
function populationsize(history)
    popsize = zeros(eltype(valtype(history)),length(collect(values(history))[1]))
    for (x,vₓ) in history
        popsize .+= vₓ
    end
    return popsize
end

function mutationload(history::Dict)
    mutload = zeros(eltype(valtype(history)),length(collect(values(history))[1]))
    for (x,vₓ) in history
        mutload .+= sum(x.genes) .* vₓ
    end
    return mutload
end

function ill_individual(history::Dict)
    ills = zeros(eltype(valtype(history)),length(collect(values(history))[1]))
    for (x,vₓ) in history
        ills .+= !x.ispropagable && vₓ
    end
    return ills
end


end  # end of module DiploidModel
