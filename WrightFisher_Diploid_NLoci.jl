#=
Implementation of a simple Wright Fisher model with constant
population size
=#

using Distributions
using SparseArrays
using Random


function WrightFisher.setupparameter(par, n0, historylength)
    #chromosome cuts of no interest for full recombination, because genes are independent in that case
    #otherwise the number of cuts is Poisson distributed, whereas the positions are uniformly choosen
    ccuts=initcuts(par)

    return (
    par...,
    rndm = Vector{Int}(undef, 2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = 1:par.Nloci,
    traits = [inittraits(par, n0),inittraits(par, n0)],
    indices = Dict(
        "healthy" => [
            collect(n0["Ill"]+1:n0["PopSize"]),
            collect(n0["Ill"]+1:n0["PopSize"]),
        ],
        "ill" => [
            collect(1:n0["Ill"]),
            collect(n0["Ill"]+1:n0["PopSize"])
        ],
    ),
    historylength = historylength,
    #chromosome cuts of no interest for full recombination, because genes are independent in that case
    #otherwise the number of cuts is Poisson distributed, whereas the positions are uniformly choosen
    ccuts = ccuts,
    choosecopy = Vector{Int64}(undef,length(ccuts)),
    choosecopyfrom = 1:2
    )
end

function WrightFisher.birth!(ps, par, indx, old, new)
    #choose two genetic configurations to mate
    rand!(par.rndm, par.indices["healthy"][old])
    for i in par.choosecopyfrom, j in par.choosecopyfrom
        dropzeros!(par.traits[old][par.rndm[i]][j])
    end
    #generate offsprings genetic configuration
    offspring!(indx, par, rand(par.MutationsPerBirth), old, new)
    #add to population size
    ps["PopSize"] += 1
    #add the individual to the current population state dictionary
    if DiploidModel2.ispropagable(par.traits[new][indx])
        #save index of new individual in appropriate list
        WrightFisher.saveorpush!(par.indices["healthy"][new], ps["PopSize"] - ps["Ill"], indx)
    else
        ps["Ill"] += 1
        #save index of new individual in approriate list
        WrightFisher.saveorpush!(par.indices["ill"][new], ps["Ill"], indx)
    end
    ps["ML"] += DiploidModel2.mutationload(par.traits[new][indx])
end

function offspring!(offspring_index, par, n_mut, old, new)
    #randomly recombine the parental genetic information
    #first for one then for the other parent
    for i in par.choosecopyfrom # =1:2
        #randomly choose one copy for each chromosome/gene block
        rand!(par.choosecopy,par.choosecopyfrom)
        for (r,chromosome) in enumerate(par.ccuts)
            view(par.traits[new][offspring_index][i],chromosome) .=
                view(par.traits[old][par.rndm[i]][par.choosecopy[r]],chromosome)
        end
    end
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ = 1:n_mut
            par.traits[new][offspring_index][rand(par.choosecopyfrom)][rand(par.MutationLocation)] = 1
        end
    end
    nothing
end

const_popsize(popsize,par,t) = par.K
random_fluctuations(popsize,par,t) = round(Integer,rand(par.fluctuations))


#---
