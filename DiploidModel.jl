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
using DataFrames
using Random

include("MainFunctions.jl")
import .Gillespie

function rungillespie(time,n₀,model_parameter)

    #setup empty population history
    l = length(time)
    population_history = Dict(x=>spzeros(valtype(n₀),l) for x in keys(n₀))

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    #execute simulation
    Gillespie.run_gillespie!(
        time,n₀,
        setupparameter(model_parameter,n₀,l),
        execute_rand!,
        rates!,
        initrates,
        population_history
        )

    return population_history

end

generatehealthypopulation(popsize,Nloci) = Dict(RandomIndividual(Nloci) => popsize)


setupparameter(par,n0,historylength) = (
    par...,
    rndm = rand(2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = DiscreteUniform(1,par.Nloci),
    parent1 = spzeros(Int8,par.Nloci),
    parent2 = spzeros(Int8,par.Nloci),
    offspring = RandomIndividual(par.Nloci),
    alltraits = collect(keys(n0)),
    popsize = [sum(values(n0)),npropagable(n0)],
    historylength = historylength
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
    if haskey(d,x)
        d[x] += diff
    else
        d[RandomIndividual(x.genes)] = diff
    end
    nothing
end

function isindictandadd!(d,x,diff,par)
    if haskey(d,x)
        d[x] += diff
    else
        indx = findfirst(y->y==x,par.alltraits)
        if !isnothing(indx)
            d[par.alltraits[indx]] = diff
        else
            push!(par.alltraits,RandomIndividual(x.genes))
            d[par.alltraits[end]] = diff
        end
    end
    nothing
end

"""
    Sum of all healty individual in the population state ps
"""
function npropagable(ps)
    null = zero(valtype(ps))
    if !isempty(ps)
        return sum(x.ispropagable ? nₓ : null for (x,nₓ) in ps)
    else
        return null
    end
end

"""
    function rates!(rates,ps,par)

        > rates is a vector of the total birth and death rate as entries
        > ps is the population state as a sparse array
        > par is a named tuple with all relevant modle parameters
"""
function truerates!(rates, ps, par)
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

function allbirthrates!(rates, ps, par)
    #calculate total population size
    #par.popsize[1] = sum(values(ps))
    #calculate size of propagable subpopulation for later use
    #par.popsize[2] = npropagable(ps)
    #linear birth for all individuals regardles of disease status
    rates[1] = par.birth * par.popsize[1]
    #uniform logistic death
    rates[2] = par.popsize[1] * par.death + par.popsize[1] * (par.popsize[1] - 1) * par.competition
    nothing
end

function choosecouple!(rndm, ps, par)
    #make it a uniform random variable in (0,total_rate)
    U = DiscreteUniform(0, par.popsize[2])
    rndm[1] = rand!(U,rndm)
    #initiate
    searchparent1 = true
    searchparent2 = true
    #choose the rate at random
    for (x, nₓ) in ps
        if x.ispropagable && nₓ > 0
            if searchparent1
                rndm[1] -= nₓ
                rndm[1] ≤ 0 && (par.parent1 .= x.genes; searchparent1 = false)
            end
            if searchparent2
                rndm[2] -= nₓ
                rndm[2] ≤ 0 && (par.parent2 .= x.genes; searchparent2 = false)
            end
        end
        !searchparent1 && !searchparent2 && break
    end
    nothing
end

function offspring!(par, n_mut)
    #randomly recombine the parental genetic information
    par.offspring.genes .= choose.(par.parent1, par.parent2)
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ in 1:n_mut
            mutate!(par,rand(par.MutationLocation))
        end
    end
    nothing
end

function mutate!(par,location)
    a = par.offspring.genes[location]
    #if there are still mutations to add
    if a < 2
        #add a mutation if a=0 and if a=1 add a mutation with prob. 1/2
        par.offspring.genes[location] = a*rand(Bool)+1
    end
end

function choose(i, j)
    if i == j == 0
        return 0
    elseif iszero(i*j)
        return rand(Bool) ? 0 : 1
    else
        return (rand(Bool) ? 1 : (rand(Bool) ? 0 : 2))
    end
end

function birth!(ps, par)
    #choose two genetic configurations to mate
    choosecouple!(par.rndm, ps, par)
    dropzeros!(par.parent1); dropzeros!(par.parent2)
    #generate offsprings genetic configuration
    offspring!(par,rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    isindictandadd!(ps,par.offspring,one(valtype(ps)))
    #increase the population size and the number of propagable individuals
    update_popsize!(par,ispropagable(par.offspring.genes),1)
end

function birth_rand!(ps,par)
    pskeys = keys(ps)
    par.parent1 .= rand(pskeys).genes
    par.parent2 .= rand(pskeys).genes
    #generate offsprings genetic configuration
    offspring!(par,rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    isindictandadd!(ps,par.offspring,one(valtype(ps)))
    #increase the population size and the number of propagable individuals
    update_popsize!(par,ispropagable(par.offspring.genes),1)
end

function update_popsize!(par,isprop,diff)
    par.popsize[1] += diff
    isprop && (par.popsize[2] += diff)
end

function death!(ps, par)
    #make it a uniform random variable in (0,pop_size)
    rndm = rand(DiscreteUniform(0,par.popsize[1]))
    #choose the rate at random and substract from population state dictionary
    for (x, nₓ) in ps
        if nₓ > 0
            rndm -= nₓ
            if rndm ≤ 0
                ps[x] -= one(valtype(ps))
                update_popsize!(par,x.ispropagable,-1)
                break
            end
        end
    end
    nothing
end

function death_rand!(ps,par)
    fey = rand(keys(ps))
    ps[fey] -= one(valtype(ps))
    update_popsize!(par,fey.ispropagable,-1)
    nothing
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

function execute_rand!(i,ps,par)
    if i == 1
        birth_rand!(ps,par)
    elseif i == 2
        death_rand!(ps,par)
    else
        error("Unknown event index: $i")
    end
end

"""
    populationsize

    Get the population size process from the history dictionary
"""
function populationsize(history)
    popsize = zeros(eltype(valtype(history)),length(last(first(history))))
    for (x,vₓ) in history
        popsize .+= vₓ
    end
    return popsize
end

function abs_mutationload(history)
    mutload = zeros(eltype(valtype(history)),length(last(first(history))))
    for (x,vₓ) in history
        mutload .+= sum(x.genes) .* vₓ
    end
    return mutload
end

function abs_ill_individual(history)
    ills = zeros(eltype(valtype(history)),length(last(first(history))))
    for (x,vₓ) in history
        ills .+= !x.ispropagable && vₓ
    end
    return ills
end

function rescale!(res_quantity,abs_quantity,abs_size)
    for (i,(x,y)) in enumerate(zip(abs_quantity,abs_size))
        res_quantity[i] = y > zero(eltype(abs_size)) ? x/y : zero(eltype(res_quantity))
    end
    nothing
end

function rescale!(abs_quantity,abs_size)
    rescale!(abs_quantity,abs_quantity,abs_size)
    nothing
end

function rescale(abs_quantity,abs_size)
    res_quantity = zeros(length(abs_quantity))
    rescale!(res_quantity,abs_quantity,abs_size)
    return res_quantity
end

function ndifferenttypes(history)
    ntypes = zeros(Int,length(last(first(history))))
    for (x,vₓ) in history
        ntypes .+= isnotzero.(vₓ)
    end
    return ntypes
end

function rel_mutation_distribution(history)
    #get any dictionary entry
    x,vₓ = first(history)
    mutations_per_time = [zeros(length(x.genes)) for _ in 1:length(vₓ)]

    for (x,vₓ) in history
    	for (time,nₓ) in zip(findnz(vₓ)...)
    		for (site,mut) in zip(findnz(x.genes)...)
    			mutations_per_time[time][site] += mut * nₓ
    		end
    	end
    end

    mutations_per_time ./= populationsize(history)

    return mutations_per_time
end

isnotzero(n) = !iszero(n)

function historytodataframe(t,Nloci,history)
    l = length(t)
    df = DataFrame(
        time=t,
        PopSize=populationsize(history),
        Ill=abs_ill_individual(history),
        Mutation=abs_mutationload(history),
        NumTypes=ndifferenttypes(history),
        )
    #generate empty mutation history
    mutationhistory = [zeros(Int,l) for _ in 1:Nloci]
    #fill mutation history
    for (x,vₓ) in history
        config = findnz(x.genes)
        timeandsize = findnz(vₓ)
        #iterate over all loci that hold mutations
        for (locus,nmut) in zip(first(config),last(config))
            #iterate over history of that type
            for (row,size) in zip(first(timeandsize),last(timeandsize))
                mutationhistory[locus][row] += nmut*size
            end
        end
    end
    #save mutation history
    for n in 1:Nloci
        df[!,"$n"] = mutationhistory[n]
    end

    return df
end

function timesfromswitches(totaltime,switches)
    times = Vector{typeof(totaltime)}(undef,length(switches)+1)
    tstart = 1
    for (i,t) in enumerate(switches)
        tend = findfirst(x -> x≥t,totaltime)
        times[i] = totaltime[tstart:tend]
        tstart = tend
    end
    times[end] = totaltime[tstart:end]
    return times
end

function runseveralgillespie(totaltime,switch::Vector,x0,parameters)
    #setup empty population history
    l = length(totaltime)
    population_history = Dict(x=>spzeros(valtype(x0),l) for x in keys(x0))
    #save initial population
    for (x,nₓ) in x0
        #add the population history for the given trait
        population_history[x][1] = nₓ
    end
    #setup empty rates vector
    initrates = Vector{typeof(parameters[1].birth)}(undef,2)
    #switching times
    times = timesfromswitches(totaltime,switch)
    #execute simulation
    for (i,t) in enumerate(times)
        Gillespie.run_gillespie!(
            t,x0,
            setupparameter(parameters[i],x0,l),
            execute!,
            rates!,
            initrates,
            population_history,
            hstart = t[1]
            )
    end
    return population_history
end

function runseveralgillespie(totaltime,switch::Number,x0,parameters)
    return runseveralgillespie(totaltime,[switch],x0,parameters)
end

end  # end of module DiploidModel
