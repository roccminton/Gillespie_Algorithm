
function DiploidModel2.setupparameter(par,n0,historylength)
    #chromosome cuts of no interest for full recombination, because genes are independent in that case
    #otherwise the number of cuts is Poisson distributed, whereas the positions are uniformly choosen
    ccuts = initcuts(par)

    return (
    par...,
    rndm = Vector{Int}(undef,2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = 1:par.Nloci,
    traits = inittraits(par,n0),
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
    choosecopyfrom = 1:2,
    )
end

function DiploidModel2.birth!(ps, par)
    #choose two genetic configurations to mate
    rand!(par.rndm,par.indices["healthy"])
    #clean up parental configurations
    for i in par.choosecopyfrom, j in par.choosecopyfrom
        dropzeros!(par.traits[par.rndm[i]][j])
    end
    #select free index for offspring
    if isempty(par.indices["free"])
        offspring_index = length(par.traits) + 1
        push!(par.traits,emptytraits(par.Nloci))
    else
        offspring_index = pop!(par.indices["free"])
    end
    #generate offsprings genetic configuration
    offspring!(offspring_index, par, rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    DiploidModel2.updateps_birth!(ps,par,offspring_index)
end

#---

"""
Sets up an empty trait vector for a Simulation with non full recobination. Therefore
2N positions are needed, for N diploid genes under consideration.
"""
function inittraits(par,n0)
    #Setup healthy genetic information
    locs = 1:par.Nloci
    #Generate healty population with some buffer for fluctuations
    traits = [emptytraits(par.Nloci) for _ in 1:round(Int,par.K + sqrt(par.K))]
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

emptytraits(Nloci,T=Bool) = [spzeros(T,Nloci),spzeros(T,Nloci)]

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

"""
    Generates the chromosome cuts depending on the recombination rate
"""
function initcuts(par)
    if par.recombination == 1
        return [i:i for i in 1:par.Nloci]
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
