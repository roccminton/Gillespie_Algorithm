
#Define Type for Population History
struct LoadPosMLP
    mlp :: Dict
    loadpos :: Dict
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadpos = Dict(
            "Healthy"=>[emptytraits(par.Nloci) for _ in 1:l],
            "Ill"=>[emptytraits(par.Nloci) for _ in 1:l]
            )
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

DiploidModel2.updatestats_death!(ps,par,index) = update_histogram!(par.cloadpos,par.traits[index],-1)
DiploidModel2.updatestats_birth!(ps,par,index) = update_histogram!(par.cloadpos,par.traits[index],+1)

#---

DiploidModel2.addstatsparameter(ph::LoadPosMLP,par,n0,l) = (par...,cloadposhist = initialhistogram(par,n0))

function initialhistogram(par,n0)
    hist = Dict("Healthy" => emptytraits(par.Nloci), "Ill" => emptytraits(par.Nloci))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_histogram!(hist,ind,1)
    end
    return hist
end

function savehistdata!(hhist,index,n0,par)
        hhist["Healthy"][index] .= par.cloadpos["Healthy"]
        hhist["Ill"][index] .= par.cloadpos["Ill"]
end

function update_histogram!(hist,ind,i)
    if DiploidModel2.ispropagable(ind)
        hist["Healthy"] .+= i .* ind
    else
        hist["Ill"] .+= i .* ind
    end
end
