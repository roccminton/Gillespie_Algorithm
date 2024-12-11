#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    loadhist :: Matrix
    loadpos :: Matrix
    snaps :: Array
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadpos = zeros(Int64,par.Nloci,l)
    loadhist = zeros(Int64,maxmutationload(par),l)
    snaps = zeros(Bool,par.Nloci,par.samplesize,length(par.savesnap))
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return MLPLoadHistLoadPos(mlp,loadhist,loadpos,snaps,par)
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
    savedata!(ph.loadpos,index,ps,par.cloadpos)
    #save the histogram History
    savedata!(ph.loadhist,index,ps,par.cloadhist)
    #safe the haploid histogramm data
    index in par.savesnap && savesnap!(ph.snaps,index,ps,par)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],-1,par.Nloci)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],+1,par.Nloci)
end

#change data for better storage
function convertforsaving(h)
    #list all the parameters worth saving
    safe_parameter=[
        "death","μ","Nloci","historylength","ccuts","recombination",
        "competition", "birth", "rates", "K", "samplesize", "savesnap"
        ]
    safe_h = Dict{String,Any}()
    #safe mlp as it is
    merge!(safe_h,h.mlp)
    #loadhist
    safe_h["LoadHist"] = h.loadhist
    #loadpos
    safe_h["LoadPos"] = h.loadpos
    #snaps
    safe_h["SnapShots"] = h.snaps
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end


#---

function DiploidModel2.addstatsparameter(par,n0,l)

    samplesize = 1000
    nsaves = 1000

    samplesize = min(samplesize,par.K)
    nsaves = min(nsaves,l)

     return (
        par...,
        samplesize = samplesize,
        sampleinds = Vector{Int64}(undef,samplesize),
        savesnap = round.(Int,range(1,l,length=nsaves)),
        cloadpos = initialloadpos(par,n0),
        cloadhist = initialloadhist(par,n0)
    )
end

function initialloadpos(par,n0)
    pos = zeros(Int64,par.Nloci)
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(pos,ind,1,par.Nloci)
    end
    return pos
end

function initialloadhist(par,n0)
    hist = zeros(Int64,maxmutationload(par))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadhist!(hist,ind,1,par.Nloci)
    end
    return hist
end

savedata!(hdata,index,n0,cdata) = (hdata[:,index] .= cdata)

function savesnap!(snaps,index,ps,par)
    t = findfirst(isequal(index),par.savesnap)
    #choose individuals at random
    rand!(par.sampleinds,vcat(par.indices["healthy"],par.indices["ill"]))
    #add snaps (we allways take the "first" gamete of any individual. Since we assume linkage equilibrium for the gametes
    #that is the same as choosing at random between the two gametes of any individual)
    for (i,ind) in enumerate(par.sampleinds)
        snaps[:,i,t] .= par.traits[ind][1]
    end
end

function update_loadpos!(pos,ind,i,Nloci)
    pos .+= i .* ind[1]
    pos .+= i .* ind[2]
end

function do_update_loadhist!(hist,load,i)
    while !checkbounds(Bool,hist,load)
        push!(hist,0)
    end
    hist[load] += i
end

function update_loadhist!(hist,ind,i,par)
    do_update_loadhist!(hist,sum(ind[1])+1,i)
    do_update_loadhist!(hist,sum(ind[2])+1,i)
end

#---

maxmutationload(Nloci,μ,K) = min(Nloci + quantile(Poisson(μ),1-1/K),2*Nloci)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = zeros(Int64,maxmutationload(par))
