using Tullio
using LoopVectorization

#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    mC :: Vector
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    mC = zeros(Float64,l)
    initialC0!(par,n₀)
    return MLPLoadHistLoadPos(mlp,mC,par)
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::MLPLoadHistLoadPos) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::MLPLoadHistLoadPos,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #calculate and save mean correlation
    index in par.savematrix && savecovmatrix!(ph.mC,index,ps,par)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1)
    updateC0!(ps,par.traits[index],-1)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1)
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
    safe_h["MeanCorr"] = h.mC
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

#---

function DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l)
    samplesize = min(1000,par.K)
    return (
        par...,
        samplesize = samplesize,
        sampleinds = Vector{Int64}(undef,samplesize),
        cloadpos = initialloadpos(par,n0),
        loadmeans = zeros(Float64,par.Nloci),
        loadstds = zeros(Float64,par.Nloci),
        corrmatrix = zeros(Float64,par.Nloci,par.Nloci),
        savematrix = round.(Int,range(1,l,length=1000)),
        )
end

function initialC0!(par,n0)
    n0["C0"] = 0
    for ind ∈ par.traits[1:n0["PopSize"]]
        updateC0!(n0,ind,1)
    end
end

function initialloadpos(par,n0)
    lp = emptytraits(par.Nloci,Int64)
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(lp,ind,1)
    end
    return lp
end

function savelistdata!(hlist,index,n0,clist)
    hlist[index] .= clist
end

function savecovmatrix!(mC,index,n0,par)
    #calculate means per position
    for i in 1:par.Nloci
        par.loadmeans[i] = sum(par.cloadpos[j][i] for j in par.choosecopyfrom)
    end
    par.loadmeans .= par.loadmeans ./ n0["PopSize"]
    #choose individuals at random
    rand!(par.sampleinds,vcat(par.indices["healthy"],par.indices["ill"]))
    #calculate std per position
    for i in 1:par.Nloci
        par.loadstds[i] = sqrt(1/par.samplesize *
            sum((sum(par.traits[ind][j][i] for j in par.choosecopyfrom) - par.loadmeans[i])^2 for ind in par.sampleinds)
            )
    end
    #calculate covariances
    for i in 1:par.Nloci-1
        for j in (i+1):par.Nloci
            par.corrmatrix[i,j] = 1/(par.loadstds[i]*par.loadstds[j]*(par.samplesize-1)) * sum(
                (sum(par.traits[ind][k][i] for k in par.choosecopyfrom) - par.loadmeans[i]) * (sum(par.traits[ind][k][j] for k in par.choosecopyfrom) - par.loadmeans[j])
            for ind in par.sampleinds)
        end
    end
    mC[index] = 2 * sum(par.corrmatrix) / (par.Nloci*(par.Nloci-1))
end

function update_loadpos!(lp,ind,i)
    lp .+= i .* ind
end

function updateC0!(n0,ind,i)
    n0["C0"] += i * iszero(sum(ind[1]))
    n0["C0"] += i * iszero(sum(ind[2]))
end


#---

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K^2)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = spzeros(Integer,maxmutationload(par))
