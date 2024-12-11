include("NonRandomMating/Toolkit.jl")

import .ToolkitNRM

using CairoMakie
using Statistics
using Distributions
using SparseArrays
using JLD
using LinearAlgebra
using DataFrames
using JSON


tend = 100_000
K = 10_000
N = 90
dni = 0.4

#Execute simulations
#---
#h = ToolkitNRM.execute_cont(K,dni,N,tend,0)

#Save simulations
#---
abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N"
#abs_path = mkpath(abs_path)
#save(abs_path * "/r=0.jld",ToolkitNRM.convertforsaving(h0))
#save(abs_path * "/r=1.jld",ToolkitNRM.convertforsaving(h1))

#Load Data
#---
d = load(abs_path * "/corrcluster1.jld")

#---
#calculate load positions and average over time
mot_loadpos(data, time, N) =
    vec(mean(data[1:N, time] .+ data[N+1:2N, time], dims = 2))

function CovToCor(M)
    iszero(M) && return M
    D = Diagonal(sqrt.(diag(M)))
    try
        DInv = inv(D)
        return DInv * M * DInv
    catch
        return M
    end
end

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)

extime(data, loadclass, birthtime, tend) = findfirst(
    iszero,
    view(data["HaploidLoadHist"], loadclass + 1, birthtime:tend),
)

function findextimes(data, tend = 0)
    iszero(tend) && (tend = data["historylength"])
    extimes = [extime(data, 0, 1, tend)]
    lc = 1
    while !isnothing(extimes[end])
        ext = extime(data, lc, extimes[end], tend)
        isnothing(ext) && return extimes
        push!(extimes, ext + extimes[end])
        lc += 1
    end
end

cluster_count(s) = iszero(length(s)) ? 0 : sum(count(!iszero, s, dims = 1))
n_cluster(s) = iszero(length(s)) ? 0 : size(s)[2]
cluster_dist(s) = iszero(length(s)) ? [] : vec(count(!iszero,s,dims=1))

#---
extimes = findextimes(d)

#start and end of plot window
ts = 14_000
te = 25_000


#---

prev = replace_NaN(d["Ill"] ./ d["PopSize"])
ml = replace_NaN(d["ML"] ./ d["PopSize"])

function get_mean_cluster(data,n0=1000)
    cluster_list = cluster_dist.(data["CorrCluster"])
    cluster = Set(cluster_list)
    times_dict = Dict(c => findall(x->x==c,cluster_list) for c in cluster)
    haploid_load_dict = Dict()
    for c in cluster
        if length(times_dict[c]) ≥ n0
            haploid_load_dict[c] = mean(data["HaploidLoadHist"][:,t] for t in times_dict[c])
        end
    end
    return haploid_load_dict
end

stringdata = JSON.json(get_mean_cluster(d))

open(abs_path * "/haploadclass_percluster.json","w") do f
    write(f,stringdata)
end
