using Tullio
using LoopVectorization

#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    initialC0!(par,n₀)
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return MLPLoadHistLoadPos(mlp,par)
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::MLPLoadHistLoadPos) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::MLPLoadHistLoadPos,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #break if C0 goes extinct
    iszero(ps["C0"]) && (ps["PopSize"] = 0)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    updateC0!(ps,par.traits[index],-1)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    updateC0!(ps,par.traits[index],+1)
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
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

#---

DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l) = par

function initialC0!(par,n0)
    n0["C0"] = 0
    for ind ∈ par.traits[1:n0["PopSize"]]
        updateC0!(n0,ind,1)
    end
end

function updateC0!(n0,ind,i)
    n0["C0"] += i * iszero(sum(ind[1]))
    n0["C0"] += i * iszero(sum(ind[2]))
end


#---
