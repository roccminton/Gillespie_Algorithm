#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    loadhist :: Dict
    loadpos :: Dict
    types :: Vector
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
    types = [emptytypes(par) for _ in 1:l]
    return MLPLoadHistLoadPos(mlp,loadhist,loadpos,types,par)
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
    #safe the types History
    savelistdata!(ph.types,index,ps,par.ctypes)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],-1,par.Nloci)
    update_types!(par.ctypes,par.traits[index],-1,par)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],+1,par.Nloci)
    update_types!(par.ctypes,par.traits[index],+1,par)
end

#change data for better storage
function convertforsaving(h)
    #list all the parameters worth saving
    safe_parameter=[
        "death","μ","Nloci","historylength","ccuts","recombination",
        "competition", "birth", "rates", "K"
        ]
    safe_h = Dict{String,Any}()
    #safe mlp as it is
    merge!(safe_h,h.mlp)
    #loadhist
    for (key,value) in h.loadhist
        safe_h["LoadHist" * key] = hcat(value...)
    end
    #loadpos
    for (key,value) in h.loadpos
        safe_h["LoadPos" * key] = hcat([vcat(d...) for d in value]...)
    end
    #types
    safe_h["Types"] = hcat(h.types...)
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

#---

function DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l)
    base2 = [
        [2^n for n in 0:par.Nloci-1],
        [2^n for n in par.Nloci:2*par.Nloci-1]]
    return (
        par...,
        cloadpos = initialloadpos(par,n0),
        cloadhist = initialloadhist(par,n0),
        ctypes = initialtypes(par,n0,base2),
        base2 = base2,
        dump = [Vector{Int64}(undef,par.Nloci),Vector{Int64}(undef,par.Nloci)]
    )
end

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

function initialtypes(par,n0,base2)
    types = emptytypes(par)
    for ind ∈ par.traits[1:n0["PopSize"]]
        types[indvtodec_nopar(ind,base2)+1] += 1
    end
    return types
end

function savehistdata!(hhist,index,n0,chist)
        hhist["Healthy"][index] .= chist["Healthy"]
        hhist["Ill"][index] .= chist["Ill"]
end

function savelistdata!(hlist,index,n0,clist)
        hlist[index] .= clist
end

function update_loadpos!(hist,ind,i,Nloci)
    if DiploidModel2.ispropagable(ind,Nloci)
        hist["Healthy"] .+= i .* ind
    else
        hist["Ill"] .+= i .* ind
    end
end

function update_loadhist!(hist,ind,i,Nloci)
    load = round(Integer,DiploidModel2.mutationload(ind)+1)
    if DiploidModel2.ispropagable(ind,Nloci)
        do_update_loadhist!(hist,"Healthy",load,i)
    else
        do_update_loadhist!(hist,"Ill",load,i)
    end
end

function do_update_loadhist!(hist,key,load,i)
    while !checkbounds(Bool, hist[key],load)
        push!(hist[key],0)
    end
    hist[key][load] += i
end

function update_types!(types,ind,i,par)
    types[indvtodec(ind,par)+1] += i
end

#---

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K^2)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = zeros(Int64,maxmutationload(par))
emptytypes(par) = zeros(Int64,ntypes(par))

ntypes(par) = 2^(2*par.Nloci)

function indvtodec(ind,par)
    for i in par.choosecopyfrom #1:2
        broadcast!(*,par.dump[i],ind[i] .* par.base2[i])
    end
    return round(Int64,sum(par.dump[1])+sum(par.dump[2]))
end
indvtodec_nopar(ind,base2) = sum(sum(i.*b for (i,b) in zip(ind,base2)))
