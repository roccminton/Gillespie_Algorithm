#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    loadhist :: Dict
    loadpos :: Dict
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadpos = Dict(
            "Healthy"=>[emptytraits(par.Nloci,Int64) for _ in 1:l],
            "Ill"=>[emptytraits(par.Nloci,Int64) for _ in 1:l]
            )
    loadhist = Dict(
            "Healthy"=>[emptyhistorgram(par) for _ in 1:l],
            "Ill"=>[emptyhistorgram(par) for _ in 1:l]
            )
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return MLPLoadHistLoadPos(mlp,loadhist,loadpos,par)
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::MLPLoadHistLoadPos) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::MLPLoadHistLoadPos,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #save the histogram History
    savehistdata!(ph.loadpos,index,ps,par.cloadpos)
    #save the histogram History
    savehistdata!(ph.loadhist,index,ps,par.cloadhist)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],-1,par.Nloci)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],+1,par.Nloci)
end

#overwrite JLD function for to save data as Matrix instead of Vec(Vec)
function convertforsaving(h)
    safe_h = Dict{String,Any}()
    #safe mlp as it is
    merge!(safe_h,h.mlp)
    #change key and converet to matrix
    for (key,value) in h.loadhist
        safe_h["LoadHist" * key] = hcat(value...)
    end
    #loadpos
    for (key,value) in h.loadpos
        safe_h["LoadPos" * key] = hcat([vcat(d...) for d in value]...)
    end
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        safe_h[String(k)] = v
    end
    return safe_h
end

#---

DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l) = (
    par...,
    cloadpos = initialloadpos(par,n0),
    cloadhist = initialloadhist(par,n0)
    )

function initialloadpos(par,n0)
    hist = Dict("Healthy" => emptytraits(par.Nloci,Int64), "Ill" => emptytraits(par.Nloci,Int64))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(hist,ind,1,par.Nloci)
    end
    return hist
end

function initialloadhist(par,n0)
    hist = Dict("Healthy" => emptyhistorgram(par), "Ill" => emptyhistorgram(par))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadhist!(hist,ind,1,par.Nloci)
    end
    return hist
end

function savehistdata!(hhist,index,n0,chist)
        hhist["Healthy"][index] .= chist["Healthy"]
        hhist["Ill"][index] .= chist["Ill"]
end

function update_loadpos!(hist,ind,i,Nloci)
    if DiploidModel2.ispropagable(ind,Nloci)
        hist["Healthy"] .+= i .* ind
    else
        hist["Ill"] .+= i .* ind
    end
end

function update_loadhist!(hist,ind,i,Nloci)
    if DiploidModel2.ispropagable(ind,Nloci)
        hist["Healthy"][round(Int64,DiploidModel2.mutationload(ind)+1)] += i
    else
        hist["Ill"][round(Int64,DiploidModel2.mutationload(ind)+1)] += i
    end
end

#---

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K^2)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = spzeros(Int64,maxmutationload(par))
