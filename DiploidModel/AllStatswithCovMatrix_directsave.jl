using SparseArrays
using LinearAlgebra

#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    haploidloadhist :: Vector
    loadpos :: Dict
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadpos = Dict(
            "Healthy"=>[spzeros(Int64,par.Nloci) for _ in 1:l],
            "Ill"=>[spzeros(Int64,par.Nloci) for _ in 1:l]
            )
    haploidloadhist = [spzeros(Int64,par.Nloci+1) for _ in 1:l]
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return MLPLoadHistLoadPos(mlp,haploidloadhist,loadpos,par)
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
    savelistdata!(ph.haploidloadhist,index,ps,par.cloadhist)
    #calculate and save covariance matrix
    index in par.savematrix && savecovmatrix(index,ps,par,ph.loadpos)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1)
    update_haploidloadhist!(par.cloadhist,par.traits[index],-1,par)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1)
    update_haploidloadhist!(par.cloadhist,par.traits[index],+1,par)
end

#change data for better storage
function convertforsaving(h)
    #list all the parameters worth saving
    safe_parameter=[
        "death","μ","Nloci","historylength","recombination",
        "competition", "birth", "K", "rate"
        ]
    safe_h = Dict{String,Any}()
    #safe mlp as it is
    merge!(safe_h,h.mlp)
    #haploid loadhist
    safe_h["HaploidLoadHist"] = h.haploidloadhist
    #loadpos
    for (key,value) in h.loadpos
        safe_h["LoadPos" * key] = value
    end
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

function save_data(abs_path,filename,h,i)
    #save floats
    save_parameter=[
        "death","μ","Nloci","historylength","recombination",
        "competition", "birth", "K"
        ]
    d_par = Dict{String,Float64}()
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ save_parameter && (d_par[key] = v)
    end
    save(abs_path *"/parameter.jld",d_par)
    #save lists
    save(abs_path * filename * "_$(i)MLP.jld",h.mlp)
    #save lists of lists
    d_llists = Dict{String,Any}()
    d_llists["HaploidLoadHist"] = h.haploidloadhist
    for (key,value) in h.loadpos
        d_llists["LoadPos" * key] = value
    end
    save(abs_path * filename * "_$i.jld",d_llists)
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
    samplesize = min(1000,par.K)
    #setup empty file for matrices
    file = jldopen(par.path * "CovMatrixs.jld","w")
    return (
        par...,
        samplesize = samplesize,
        samplegamets = Vector{Int64}(undef,samplesize),
        sampleinds = Vector{Int64}(undef,samplesize),
        cloadpos = initialloadpos(par,n0),
        cloadhist = initialhaploidloadhist(par,n0),
        loadmeans = zeros(Float64,par.Nloci),
        savematrix = round.(Int,range(1,l,length=1000)),
        covmatrix = zeros(par.Nloci,par.Nloci),
        file = file
        )
end

function initialloadpos(par,n0)
    hist = Dict("Healthy" => spzeros(Int64,par.Nloci), "Ill" => spzeros(Int64,par.Nloci))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(hist,ind,1)
    end
    return hist
end

function initialhaploidloadhist(par,n0)
    hist = zeros(par.Nloci+1)
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_haploidloadhist!(hist,ind,1,par)
    end
    return hist
end

function savehistdata!(hhist,index,n0,chist)
        hhist["Healthy"][index] .= chist["Healthy"]
        hhist["Ill"][index] .= chist["Ill"]
end

function savelistdata!(hlist,index,n0,clist)
    while length(hlist[index]) < length(clist)
        push!(hlist[index],zero(eltype(hlist[index])))
    end
    hlist[index] .= clist
end

function savecovmatrix(index,n0,par,loadpos)
    #calculate means per position
    for i in 1:par.Nloci
        par.loadmeans[i] = loadpos["Healthy"][index][i] + loadpos["Ill"][index][i]
    end
    par.loadmeans .= par.loadmeans ./ n0["PopSize"]
    #choose individuals at random
    rand!(par.sampleinds,vcat(par.indices["healthy"],par.indices["ill"]))
    #choose gamets at random
    rand!(par.samplegamets,par.choosecopyfrom)
    #calculate covariances
    par.covmatrix .= 0.0
    for (ind_index,i_choice) in zip(par.sampleinds,par.samplegamets)
        par.covmatrix .+= AAt!(par.covmatrix,par.traits[ind_index][i_choice].-par.loadmeans)
    end
    #C = sum(AA′(par.traits[ind_index][i_choice].-par.loadmeans) for (ind_index,i_choice) in zip(par.sampleinds,par.samplegamets))
    #normalize
    par.covmatrix .= par.covmatrix ./ (n0["PopSize"]-1)
    #safe
    write(par.file, "$index", par.covmatrix)
end

# calculation of A*A' for a Matrix A
AA′(A) = A*A'
AAt!(C,A) = mul!(C,A,transpose(A))

function update_loadpos!(hist,ind,i)
    if DiploidModel2.ispropagable(ind)
        hist["Healthy"] .+= i .* sum(ind)
    else
        hist["Ill"] .+= i .* sum(ind)
    end
end

function update_haploidloadhist!(hist,ind,i,par)
    for j in par.choosecopyfrom #1:2
        hist[sum(ind[j])+1] += i
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
