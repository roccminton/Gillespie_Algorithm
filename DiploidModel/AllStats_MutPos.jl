"""Saves MLP, Load Histogram and Mutation Positions split into double and single
mutations. Can only be used together with "BirthDeathwithRec" because we need
the two copy structure. Could be modified for other structures as well."""

#Define Type for Population History
struct MLPLoadHistMutPos
    mlp :: Dict
    loadhist :: Dict
    mutpos :: Dict
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    mutpos = Dict(
            "Single"=>[spzeros(par.Nloci) for _ in 1:l],
            "Double"=>[spzeros(par.Nloci) for _ in 1:l]
            )
    loadhist = Dict(
            "Healthy"=>[emptyhistorgram(par) for _ in 1:l],
            "Ill"=>[emptyhistorgram(par) for _ in 1:l]
            )
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return MLPLoadHistMutPos(mlp,loadhist,mutpos,par)
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::MLPLoadHistMutPos) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::MLPLoadHistMutPos,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #save the histogram History
    savehistdata!(ph.mutpos,index,ps,par.cmutpos)
    #save the histogram History
    savehistdata!(ph.loadhist,index,ps,par.cloadhist)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_mutpos!(par.cmutpos,par.traits[index],-1)
    update_loadhist!(par.cloadhist,par.traits[index],-1)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_mutpos!(par.cmutpos,par.traits[index],+1)
    update_loadhist!(par.cloadhist,par.traits[index],+1)
end

#---

DiploidModel2.addstatsparameter(ph::MLPLoadHistMutPos,par,n0,l) = (
    par...,
    cmutpos = initialmutpos(par,n0),
    cloadhist = initialloadhist(par,n0)
    )

function initialmutpos(par,n0)
    hist = Dict("Single" => spzeros(par.Nloci), "Double" => spzeros(par.Nloci))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_mutpos!(hist,ind,1)
    end
    return hist
end

function initialloadhist(par,n0)
    hist = Dict("Healthy" => emptyhistorgram(par), "Ill" => emptyhistorgram(par))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadhist!(hist,ind,1)
    end
    return hist
end

function savehistdata!(hhist,index,n0,chist)
    for key in keys(hhist)
        hhist[key][index] .= chist[key]
    end
end


function update_mutpos!(hist,ind,i)
    hist["Single"] .+= i .* abs.(ind[1] .- ind[2])
    hist["Double"] .+= i .* (ind[1] .* ind[2])
end

function update_loadhist!(hist,ind,i)
    if DiploidModel2.ispropagable(ind)
        hist["Healthy"][round(Integer,DiploidModel2.mutationload(ind)+1)] += i
    else
        hist["Ill"][round(Integer,DiploidModel2.mutationload(ind)+1)] += i
    end
end

#---

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = spzeros(Integer,maxmutationload(par))
