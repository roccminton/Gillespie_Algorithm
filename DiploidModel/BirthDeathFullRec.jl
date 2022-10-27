
function DiploidModel2.setupparameter(par,n0,historylength)

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
    )
end

function DiploidModel2.birth!(ps, par)
    #choose two genetic configurations to mate
    rand!(par.rndm,par.indices["healthy"])
    dropzeros!(par.traits[par.rndm[1]]); dropzeros!(par.traits[par.rndm[2]])
    #select free index for offspring
    if isempty(par.indices["free"])
        offspring_index = length(par.traits) + 1
        push!(par.traits,emptytrait(par.Nloci))
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
Sets up an empty trait vector for a Simulation with independent genes (full recobination)
Only N positions are needed, if N diploid (!) genes are considered.
"""
function inittraits(par,n0)
    locs = 1:par.Nloci
    traits = [emptytrait(par.Nloci) for _ in 1:round(Int,par.K + sqrt(par.K))]
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

emptytrait(Nloci) = spzeros(Nloci)

function offspring!(offspring_index, par, n_mut)
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
