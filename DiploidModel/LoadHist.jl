
#Define Type for Population History
struct LoadHistMLP
    mlp :: Dict
    loadhist :: Dict
    par :: NamedTuple
end

#Setup an empty population history to fill during the runtime of the simulation
function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadhist = Dict(
            "Healthy"=>[emptyhistorgram(par) for _ in 1:l],
            "Ill"=>[emptyhistorgram(par) for _ in 1:l]
            )
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return LoadHistMLP(mlp,loadhist,par)
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::LoadHistMLP) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::LoadHistMLP,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #save the histogram History
    savehistdata!(ph.loadhist,index,ps,par)
end

DiploidModel2.updatestats_death!(ps,par,index) = update_histogram!(par.cloadhist,par.traits[index],-1)
DiploidModel2.updatestats_birth!(ps,par,index) = update_histogram!(par.cloadhist,par.traits[index],+1)

#---

DiploidModel2.addstatsparameter(ph::LoadHistMLP,par,n0,l) = (par...,cloadhist = initialhistogram(par,n0))

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = spzeros(Integer,maxmutationload(par))

function initialhistogram(par,n0)
    hist = Dict("Healthy" => emptyhistorgram(par), "Ill" => emptyhistorgram(par))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_histogram!(hist,ind,1)
    end
    return hist
end

function savehistdata!(hhist,index,n0,par)
        hhist["Healthy"][index] .= par.cloadhist["Healthy"]
        hhist["Ill"][index] .= par.cloadhist["Ill"]
end

function update_histogram!(hist,ind,i)
    if DiploidModel2.ispropagable(ind)
        hist["Healthy"][round(Integer,DiploidModel2.mutationload(ind)+1)] += i
    else
        hist["Ill"][round(Integer,DiploidModel2.mutationload(ind)+1)] += i
    end
end
