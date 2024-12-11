using Tullio
using LoopVectorization

#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    loadhist :: Dict
    loadpos :: Dict
    hamdistvec :: Vector
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
    hamdistvec = Vector{Matrix{Integer}}(undef,l)
    return MLPLoadHistLoadPos(mlp,loadhist,loadpos,hamdistvec,par)
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
    #calculate and save covariance matrix
    savehammingdistance!(ph.hamdistvec,index,ps,par)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1)
    update_loadhist!(par.cloadhist,par.traits[index],-1)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1)
    update_loadhist!(par.cloadhist,par.traits[index],+1)
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
        levelhists!(value)
        safe_h["LoadHist" * key] = hcat(value...)
    end
    #loadpos
    for (key,value) in h.loadpos
        safe_h["LoadPos" * key] = hcat([vcat(d...) for d in value]...)
    end
    #haploid loadhist
    safe_h["HammDistList"] = h.hamdistvec
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

function levelhists!(histlist)
    maxlength = maximum(length.(histlist))
    for hist in histlist
        while length(hist)<maxlength
            push!(hist,zero(eltype(hist)))
        end
    end
end
#---

function DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l)
    samplesize = 1000
    return (
        par...,
        samplesize = samplesize,
        sampleindices = Vector{Integer}(undef,samplesize),
        samplechromosome = rand(par.choosecopyfrom,samplesize),
        cloadpos = initialloadpos(par,n0),
        cloadhist = initialloadhist(par,n0),
        loadmeans = zeros(Float64,par.Nloci),
        )
end

function initialloadpos(par,n0)
    hist = Dict("Healthy" => emptytraits(par.Nloci,Int64), "Ill" => emptytraits(par.Nloci,Int64))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(hist,ind,1)
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
        hhist["Healthy"][index] .= chist["Healthy"]
        hhist["Ill"][index] .= chist["Ill"]
end

function savehammingdistance!(hamdistvec,index,n0,par)
    #initialize hamming distance matrix
    H = Matrix{Integer}(undef,(2,par.samplesize))
    #choose random individual for comparison
    rand!(par.sampleindices,vcat(par.indices["healthy"],par.indices["ill"]))
    #set individual to compare with
    compareind_index = par.sampleindices[1]
    #iterate through individuals
    for (n,(ind_index,chroms)) in enumerate(zip(par.sampleindices,par.samplechromosome))
        for j in par.choosecopyfrom
            H[j,n] = hammingdistance(ind_index,compareind_index,chroms,j,par)
        end
    end
    #safe
    hamdistvec[index] = H
end

#calculation of hamming distance
hammingdistance(ind_index,compareind_index,i,j,par) = sum(x != y for (x,y) in zip(par.traits[ind_index][i],par.traits[compareind_index][j]))


function update_loadpos!(hist,ind,i)
    if DiploidModel2.ispropagable(ind)
        hist["Healthy"] .+= i .* ind
    else
        hist["Ill"] .+= i .* ind
    end
end

function update_loadhist!(hist,ind,i)
    load = round(Integer,DiploidModel2.mutationload(ind)+1)
    if DiploidModel2.ispropagable(ind)
        update_loadhist!(hist,"Healthy",load,i)
    else
        update_loadhist!(hist,"Ill",load,i)
    end
end

function update_loadhist!(hist,key,load,i)
    while !checkbounds(Bool, hist[key],load)
        push!(hist[key],0)
    end
    hist[key][load] += i
end

#---

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K^2)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = spzeros(Integer,maxmutationload(par))
