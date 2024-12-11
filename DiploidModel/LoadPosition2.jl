
#Define Type for Population History
struct LoadPosMLP
    mlp :: Dict
    loadpos :: Array
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    #3D Matrix where we have the following dimensions (#Mutations per Locus (0,1,2), Locus (1-N), time)
    loadpos = zeros(Integer,(3,par.Nloci,l))
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return LoadPosMLP(mlp,loadpos,par)
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::LoadPosMLP) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::LoadPosMLP,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #save the histogram History
    savehistdata!(ph.loadpos,index,ps,par)
end

DiploidModel2.updatestats_death!(ps,par,index) = update_histogram!(par.cloadposhist,par.traits[index],-1,par)
DiploidModel2.updatestats_birth!(ps,par,index) = update_histogram!(par.cloadposhist,par.traits[index],+1,par)

#---

DiploidModel2.addstatsparameter(ph::LoadPosMLP,par,n0,l) = (
    par...,
    cloadposhist = initialhistogram(par,n0)
    )

function initialhistogram(par,n0)
    hist = zeros(Integer,(3,par.Nloci))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_histogram!(hist,ind,1,par)
    end
    return hist
end

function savehistdata!(hhist,index,n0,par)
        view(hhist,:,:,index) .= par.cloadposhist
end

function update_histogram!(hist,ind,i,par)
    for n in 1:par.Nloci
        hist[diploid_mutation(ind,n)+1,n] += i
    end
end

diploid_mutation(ind::Vector,n) = ind[1][n] + ind[2][n]
diploid_mutation(ind::SparseVector,n) = ind[n]
