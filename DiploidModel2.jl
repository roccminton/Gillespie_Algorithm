"""

module DiploidModel2

Implementation of the same DiploidModel, but the genetic information from the
individuals in the population is saved only once in shared memory and then only
overwritten. Therefore 3 Arrays are implemented containing the indexes of the
individuals that are 1) alive and healthy 2) alive and ill and 3) not alive, hence
free space to use

birth = b,
death = d,
competition = (b-d)/K,
μ = dni,
Nloci = N,
K = K,
rates = "allbirthrates!"
"""

module DiploidModel2

using SparseArrays
using Random
using Distributions
using DataFrames
using StatsBase

include("MainFunctions.jl")
import .Gillespie


"""
    Executes one run of the Gillespie algorithm with initial population size n0 over
    the time intervall time. The named tuple model_parameter has to have at least the following
    keys: birth, death, competition, μ, Nloci, K, recombination
"""
function rungillespie(time,n₀,model_parameter)

    l = length(time)

    #setup empty population history
    #if ther is an inital histogram in the parameters than use a different history
    if haskey(model_parameter, :histogramdata)
        population_history = Dict{String,Any}(
            "Healthy"=>[spzeros(Integer,maxmutationload(model_parameter)) for _ in 1:l],
            "Ill"=>[spzeros(Integer,maxmutationload(model_parameter)) for _ in 1:l]
            )
    else
        population_history = Dict{String,Any}(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    end

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel2,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    parameter = setupparameter(model_parameter,n₀,l)

    #execute simulation
    if haskey(model_parameter, :histogramdata)
        Gillespie.run_gillespie!(
            time,n₀,
            parameter,
            (model_parameter.recombination == 1 ? execute_fullrec! : execute!),
            rates!,
            initrates,
            population_history,
            statistic! = savehistdata!
            )
    else
        Gillespie.run_gillespie!(
            time,n₀,
            parameter,
            (model_parameter.recombination == 1 ? execute_fullrec! : execute!),
            rates!,
            initrates,
            population_history,
            )
    end

    #save additional information
    population_history["cutsat"] = vcat([1],[cut[end] for cut in parameter.ccuts])

    return population_history
end

generatehealthypopulation(popsize) = initpopulation(popsize,0,0)

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

function initpopulation(popsize,ML,Ill)
    if ML >= 2*Ill
        return Dict(
            "PopSize" => popsize,
            "Ill" => Ill,
            "ML" => ML
            )
    else
        error("Mutation Load ($ML) must be at least twice the number of ill individual($Ill)")
    end
end

function initialhistogram(par,n₀)
    for ind ∈ 1:n₀["Ill"]
        par.histogramdata["Ill"][mutationload(par.traits[ind])] += 1
    end
    for ind ∈ n₀["Ill"]+1:n₀["PopSize"]
        par.histogramdata["Healthy"][mutationload(par.traits[ind])] += 1
    end
end

emptyhistogram(Nloci,μ,K) = Dict(
    "Healthy" => spzeros(Integer,maxmutationload(Nloci,μ,K)),
    "Ill" => spzeros(Integer,maxmutationload(Nloci,μ,K)),
)

function healthyhistogram(Nloci,μ,K)
    hist = emptyhistogram(Nloci,μ,K)
    hist["Healthy"][1] = K
    return hist
end

#---
#General Model Functions for parameter changes during runtime

function changeparameter(par,newpar)
    return merge(par,newpar)
end

function runseveralgillespie(totaltime,switch::Vector,parameter_changes,n₀,model_parameter)

    #setup empty population history
    l = length(totaltime)
    population_history = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel2,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    par = setupparameter(model_parameter,n₀,l)

    #save initial population
    for (x,nₓ) in n₀
        #add the population history for the given trait
        population_history[x][1] = nₓ
    end

    #switching times
    times = timesfromswitches(totaltime,switch)

    #execute simulation
    for (i,t) in enumerate(times)
        Gillespie.run_gillespie!(
            t,n₀,
            par,
            execute!,
            rates!,
            initrates,
            population_history,
            hstart = t[1]
            )
        #change parameter for next block, unless no block follows
        i < length(times) && (par = changeparameter(par,parameter_changes[i]))
    end
    return population_history
end

function runseveralgillespie(totaltime,switch::Number,par_changes,x0,parameters)
    return runseveralgillespie(totaltime,[switch],par_changes,x0,parameters)
end

#---
#Model Functions - Total Recombination

function setupparameter(par,n0,historylength)
    #chromosome cuts of no interest for full recombination, because genes are independent in that case
    #otherwise the number of cuts is Poisson distributed, whereas the positions are uniformly choosen
    ccuts = initcuts(par)

    return (
    par...,
    rndm = Vector{Int}(undef,2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = 1:par.Nloci,
    traits = inittraitsfromdict(par,n0),
    indices = Dict(
        "healthy" => collect(n0["Ill"]+1:n0["PopSize"]),
        "ill" => collect(1:n0["Ill"]),
        "free" => collect(n0["PopSize"]+1:round(Int,par.K + sqrt(par.K)))
    ),
    historylength = historylength,
    #chromosome cuts of no interest for full recombination, because genes are independent in that case
    #otherwise the number of cuts is Poisson distributed, whereas the positions are uniformly choosen
    ccuts = ccuts,
    choosecopy = Vector{Int64}(undef,length(ccuts)),
    choosecopyfrom = 1:2
    )
end

"""
    Generates the chromosome cuts depending on the recombination rate
"""
function initcuts(par)
    if par.recombination == 1
        return []
    else
        cutsat = sort!(sample(1:par.Nloci-1,rand(Poisson(par.recombination*par.Nloci)),replace=false))
        isempty(cutsat) && return [1:par.Nloci]
        ccuts = [1:cutsat[1]]
        for i in 2:length(cutsat)
            push!(ccuts,cutsat[i-1]+1:cutsat[i])
        end
        push!(ccuts,cutsat[end]+1:par.Nloci)
        return ccuts
    end
end
"""
    Choose the respective initial traits depending on the recombination rate.
"""
function inittraitsfromdict(par,n0)
    if par.recombination == 1
        return inittraits_fullrec(par,n0)
    else
        return inittraits_rec(par,n0)
    end
end

"""
Sets up an empty trait vector for a Simulation with independent genes (full recobination)
Only N positions are needed, if N diploid (!) genes are considered.
"""
function inittraits_fullrec(par,n0)
    locs = 1:par.Nloci
    traits = [spzeros(par.Nloci) for _ in 1:round(Int,par.K + sqrt(par.K))]
    for i in 1:n0["Ill"]
        l = rand(locs)
        traits[i][l] = 2
    end
    individuals = n0["Ill"]+1:n0["PopSize"]
    for i in n0["Ill"]+1:n0["ML"]-2*n0["Ill"]
        ind = rand(individuals)
        l = rand(locs)
        while traits[ind][l] ≠ 0
            ind = rand(individuals)
            l = rand(locs)
        end
        traits[ind][l] = 1
    end
    return traits
end
"""
Sets up an empty trait vector for a Simulation with non full recobination. Therefore
2N positions are needed, for N diploid genes under consideration.
"""
function inittraits_rec(par,n0)
    #Setup healthy genetic information
    locs = 1:par.Nloci
    #Generate healty population with some buffer for fluctuations
    traits = [(spzeros(par.Nloci),spzeros(par.Nloci)) for _ in 1:round(Int,par.K + sqrt(par.K))]
    #add two mutations to completely healthy individuals to get the required number of ill individuals
    for i in 1:n0["Ill"]
        l = rand(locs)
        traits[i][1][l] = 1
        traits[i][2][l] = 1
    end
    individuals = 1:n0["PopSize"]
    #add the remaining mutaions to the population to get the required mutation load
    for i in n0["Ill"]+1:n0["ML"]-2*n0["Ill"]
        #choose random individual and location
        ind = rand(individuals)
        l = rand(locs)
        #recoose random individual and location if the individual has already a mutation
        #at that locationo or at the homologe gene
        while traits[ind][1][l]+traits[ind][2][l] ≠ 0
            ind = rand(individuals)
            l = rand(locs)
        end
        traits[ind][rand(par.choosecopyfrom)][l] = 1
    end
    return traits
end

"""
    function rates!(rates,ps,par)

        > rates is a vector of the total birth and death rate as entries
        > ps is the population state as a sparse array
        > par is a named tuple with all relevant modle parameters
"""
function truerates!(rates, ps, par)
    #linear birth for all propagable individuals
    rates[1] = par.birth * (ps["PopSize"]-ps["Ill"])
    #uniform logistic death
    rates[2] = ps["PopSize"] * par.death + ps["PopSize"] * (ps["PopSize"] - 1) * par.competition
    nothing
end

function allbirthrates!(rates, ps, par)
    #linear birth for all individuals
    rates[1] = (par.birth * ps["PopSize"]) * (!isempty(par.indices["healthy"]))
    #uniform logistic death
    rates[2] = ps["PopSize"] * par.death + ps["PopSize"] * (ps["PopSize"] - 1) * par.competition
    nothing
end

#--- Functions for the model with 100% Recombination (full recombination)

function birth_fullrec!(ps, par)
    #choose two genetic configurations to mate
    rand!(par.rndm,par.indices["healthy"])
    dropzeros!(par.traits[par.rndm[1]]); dropzeros!(par.traits[par.rndm[2]])
    #select free index for offspring
    if isempty(par.indices["free"])
        offspring_index = length(par.traits) + 1
        push!(par.traits,spzeros(par.Nloci))
    else
        offspring_index = pop!(par.indices["free"])
    end
    #generate offsprings genetic configuration
    offspring_fullrec!(offspring_index, par, rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    updateps_birth!(ps,par,offspring_index)
    #update further statistics
    par.updatestats_birth!(ps,par,offspring_index)
end

function offspring_fullrec!(offspring_index, par, n_mut)
    #randomly recombine the parental genetic information
    par.traits[offspring_index] .= choose.(
        par.traits[par.rndm[1]],
        par.traits[par.rndm[2]]
    )
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ in 1:n_mut
            mutate!(par,offspring_index,rand(par.MutationLocation))
        end
    end
    nothing
end

function mutate!(par,offspring_index,location)
    a = par.traits[offspring_index][location]
    #if there are still mutations to add
    if a < 2
        #add a mutation if a=0 and if a=1 add a mutation with prob. 1/2
        par.traits[offspring_index][location] = executemutation(a)
    end
end

function executemutation(existing_mutations)
    if iszero(existing_mutations)
        return 1
    else
        rand(Bool) ? (return 1) : (return 2)
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

ispropagable(a::SparseVector) = !(2 ∈ a.nzval)

function execute_fullrec!(i,ps,par)
    if i == 1
        birth_fullrec!(ps,par)
    elseif i == 2
        death!(ps,par)
    else
        error("Unknown event index: $i")
    end
end


##--- Functions for the model with moderate to non recombination
function birth!(ps, par)
    #choose two genetic configurations to mate
    rand!(par.rndm,par.indices["healthy"])
    #clean up parental configurations
    for i in par.choosecopyfrom, j in par.choosecopyfrom
        dropzeros!(par.traits[par.rndm[i]][j])
    end
    #select free index for offspring
    if isempty(par.indices["free"])
        offspring_index = length(par.traits) + 1
        push!(par.traits,(spzeros(par.Nloci),spzeros(par.Nloci)))
    else
        offspring_index = pop!(par.indices["free"])
    end
    #generate offsprings genetic configuration
    offspring!(offspring_index, par, rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    updateps_birth!(ps,par,offspring_index)
    #update further statistics
    par.updatestats_birth!(ps,par,offspring_index)
end

function updateps_birth!(ps,par,offspring_index)
    if ispropagable(par.traits[offspring_index])
        push!(par.indices["healthy"],offspring_index)
    else
        ps["Ill"] += 1
        push!(par.indices["ill"],offspring_index)
    end
    ps["PopSize"] += 1
    ps["ML"] += mutationload(par.traits[offspring_index])
end

function updateps_death(ps,par,fey_index)
    ps["PopSize"] -= 1
    ps["ML"] -= mutationload(par.traits[fey_index])
    nothing
end

function savehistdata!(pop_hist,index,n0,par)
        pop_hist["Healthy"][index] .= par.histogramdata["Healthy"]
        pop_hist["Ill"][index] .= par.histogramdata["Ill"]
end

nostats(ps,par) = nothing

function update_histogram!(ps,par,index,i)
    if ispropagable(par.traits[index])
        par.histogramdata["Healthy"][round(Integer,mutationload(par.traits[index])+1)] += i
    else
        par.histogramdata["Ill"][round(Integer,mutationload(par.traits[index])+1)] += i
    end
end

update_histogram_birth!(ps,par,index) = update_histogram!(ps,par,index,+1)
update_histogram_death!(ps,par,index) = update_histogram!(ps,par,index,-1)

function offspring!(offspring_index, par, n_mut)
    #randomly recombine the parental genetic information
    #first for one then for the other parent
    for i in par.choosecopyfrom # =1:2
        #randomly choose one copy for each chromosome/gene block
        rand!(par.choosecopy,par.choosecopyfrom)
        for (r,chromosome) in enumerate(par.ccuts)
            view(par.traits[offspring_index][i],chromosome) .=
                view(par.traits[par.rndm[i]][par.choosecopy[r]],chromosome)
        end
    end
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ in 1:n_mut
            par.traits[offspring_index][rand(par.choosecopyfrom)][rand(par.MutationLocation)] = 1
        end
    end
    nothing
end

function ispropagable(a::Tuple{SparseVector,SparseVector})
    for (i,gene) in enumerate(a[1])
        isone(gene) && isone(a[2][i]) && return false
    end
    return true
end

mutationload(a::Tuple{SparseVector,SparseVector}) = sum(sum(spvec) for spvec in a)
mutationload(a::SparseVector) = sum(a)

function death!(ps, par)
    #choose fey
    if rand()<=ps["Ill"]/ps["PopSize"]
        fey_index = popat!(par.indices["ill"],rand(1:ps["Ill"]))
        ps["Ill"] -= 1
    else
        fey_index = popat!(par.indices["healthy"],rand(1:(ps["PopSize"]-ps["Ill"])))
    end
    #add fey to gravejard
    push!(par.indices["free"],fey_index)
    #update population state
    updateps_death(ps,par,fey_index)
    #update further statistics
    par.updatestats_death!(ps,par,fey_index)
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

isnotzero(n) = !iszero(n)

historytodataframe(history) = DataFrame(
    PopSize=history["PopSize"],
    Ill=history["Ill"],
    Mutation=history["ML"],
    Ccuts = vcat(history["cutsat"],fill(NaN,length(history["PopSize"])-length(history["cutsat"])))
    )

function timesfromswitches(totaltime,switches)
    times = Vector{typeof(totaltime)}(undef,length(switches)+1)
    tstart = 1
    for (i,t) in enumerate(switches)
        tend = findfirst(x -> x≥t,totaltime)
        times[i] = totaltime[tstart:tend]
        tstart = tend+1
    end
    times[end] = totaltime[tstart:end]
    return times
end



end  # end of module DiploidModel2
