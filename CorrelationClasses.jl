include("NonRandomMating/Toolkit.jl")

import .ToolkitNRM

using CairoMakie
using Statistics
using Distributions
using SparseArrays
using JLD
using LinearAlgebra
using ProgressMeter

tend = 100_000
K = 10_000
N = 60
dni = 0.6

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
#d = load(abs_path * "/corrcluster1.jld")

#---
#calculate load positions and average over time
mot_loadpos(data, time, N) =
    vec(mean(data[1:N, time] .+ data[N+1:2N, time], dims = 2))

function CovToCor(M)
    iszero(M) && return M
    D = Diagonal(sqrt.(diag(M)))
    DInv = inv(D)
    return DInv * M * DInv
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

function findcorrelationcluster(corrmx,ε=0.1)
    N = size(corrmx)[1]
    C = []
    for i in 1:N-1
        for j in i+1:N-1
            if corrmx[i,j] ≥ 1 - ε
                ci = findindinset(C,i)
                iszero(ci) ? push!(C,[i,j]) : (!(j ∈ C[ci]) && push!(C[ci],j))
            end
        end
    end
    return C
end

#Look in a List of Sets for an integer i.
#Returns the index of the resprective set in the list if i is in the set
#and zero if i couldnt be found
function findindinset(S,i)
    for (n,s) in enumerate(S)
        i ∈ s && return n
    end
    return 0
end

function get_scatter_data(d)
    T = [[Int[] for _ in 1:maxclust] for _ in 1:N]
    for t in 1:d["historylength"]
        for (c,col) in enumerate(eachcol(d["CorrCluster"][t]))
            for n in col
                !iszero(n) && push!(T[n][c],t)
            end
        end
    end
    return T
end

Load(N,μ,s=0,k=1) = 2N*sqrt(1-exp(-μ/N))
Prev(μ) = 1-exp(-μ)
h(N,μ) = sqrt(1-exp(-μ/N))

function Prev(N,μ,s::Matrix)
    k = size(s)[2]
    n = count(!iszero,s,dims=1)

    P=0.0
    h=sqrt(1-exp(-μ/N))
    sn = sum(n)

    for i ∈ 1:k
        for j ∈ (i+1):k
            P += (1-h)^(n[i]+n[j]) * (1-h^2)^(N-n[i]-n[j])
        end
    end

    return 1-(1-1/k^2)*(2/(k*(k-1)))*P
end

function Load(N,μ,s::Matrix)
    k = size(s)[2]
    n = count(!iszero,s,dims=1)

    L=0.0
    h=sqrt(1-exp(-μ/N))
    sn = sum(n)

    for i ∈ 1:k
        for j ∈ (i+1):k
            L += n[i]+n[j]+(2N-sn)*sqrt(1-exp(-μ/N))
        end
    end

    return (2L)/(k*(k-1))
end

#indicates when some cluster changes
function clusterchangetimes(d,f)
    T = Int64[]
    a = f(d["CorrCluster"][1])
    for (t,s) in enumerate(d["CorrCluster"])
        b = f(s)
        a ≠ b && push!(T,t)
        a = b
    end
    return T
end





#---
extimes = findextimes(d)
ctimes = clusterchangetimes(d,cluster_count)

#start and end of plot window
ts = 0
te = tend

t0,t1,t2 = 75_000, 80_000, 100_000
times = [t0:t1,t1:t2]

t = 90_000

#---

prev = replace_NaN(d["Ill"] ./ d["PopSize"])
ml = replace_NaN(d["ML"] ./ d["PopSize"])

cormatrixs = [
    CovToCor(mean(d["CovMatrixList"][int])) for int in times
]


#---
#create figure

f = Figure()

axprev = Axis(
    f[1:2, 1:4],
    yticklabelcolor = :orange,
    yaxisposition = :right,
    xaxisposition = :top,
    xticksvisible = false,
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    ylabel = "Prevalence",
    xticks = (vcat(extimes,[t]),vcat(fill(L"\dagger",length(extimes)),[L"t^\star"]))
)

axload = Axis(
    f[1:2, 1:4],
    yticklabelcolor = :red,
    ylabel = "Mutation Load",
    #xlabel = "Time",
)

axload.xticks = (
    range(ts, te; length = 5),
    [string.(round.(Integer, x / 1000)) .* "K" for x in range(ts,te;length=5)]
    )

xlims!(axload,(ts,te))
xlims!(axprev,(ts,te))
ylims!(axload,low=0)
ylims!(axprev,low=0)
linkxaxes!(axload,axprev)

#plot the data
prevline = lines!(
    axprev,
    ts:te,
    prev[ts+1:te+1],
    color = :orange,
    label = "Prevalence",
)
loadline = lines!(
    axload,
    ts:te,
    ml[ts+1:te+1],
    color = :red,
    label = "Mutation Load",
)

# hlines!(axload,[
#     Load(N,dni/2),
#     Load(N,dni/2,d["CorrCluster"][t])
#     ])
vlines!(axload,[t],color=:gray,linestyle=:dash)
#hlines!(axprev,Prev(dni/2))



#---

loadpos = mean(d["LoadPos2"][:,:,t:end],dims=3) ./ K

ax2 = Axis(f[3,1:4],title=L"Average over $ t \in [t^\star,t_{\mathrm{end}}]$")
barplot!(ax2,loadpos[1,:] .+ loadpos[2,:] .+ loadpos[3,:],color =(:gray,1))#,label="0/0")
barplot!(ax2,loadpos[2,:] .+ loadpos[3,:],color =:darkblue)#,label="0/1")
barplot!(ax2,loadpos[3,:],color =:salmon)#,label="1/1")

# ylims!(ax2,0.41,0.48); vfill = 0.475
ylims!(ax2,low=0); vfill = 1.1
xlims!(ax2,0.5,N+0.5)

color = [:red,:green,:blue,:purple,:yellow]

for i in 1:size(d["CorrCluster"][t])[2]
    c = filter(!iszero,d["CorrCluster"][t][:,i])
    scatter!(c,fill(vfill,length(c)),color=color[i],label="$(length(c))")
end

#hlines!(ax2,2*h(N,dni/2)+h(N,dni/2)^2,color=:orange)

Legend(f[4,4],ax2,framevisible=false,L"\text{# Genes}")

#---

ccill = filter(!iszero,mean(d["ClusterCountIll"][:,t:end],dims=2)) ./ 2K
cchealthy = filter(!iszero,mean(d["ClusterCountHealthy"][:,t:end],dims=2)) ./ 2K

axh = Axis(f[4,1],xlabel = "Healthy Individual",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%"
    )
axi = Axis(f[4,2],xlabel = "Ill Individual",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%"
    )
axg = Axis(f[4,3],xlabel = "Total",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%"
    )

barplot!(axi,ccill,color=[:red,:green,:blue,:purple])
barplot!(axh,cchealthy,color=[:red,:green,:blue,:purple])
barplot!(axg,ccill .+ cchealthy,color=[:red,:green,:blue,:purple])

linkyaxes!(axh,axi,axg)
for ax in [axh,axi,axg]
    ylims!(ax,low=0)
    hidexdecorations!(ax,label=false)
    ax != axh && hideydecorations!(ax,grid=false)
end

# ylims!(axh,0.26,0.3)
# ylims!(axi,0.195,0.23)
# ylims!(axg,0.46,0.53)


#---

T = round.(Integer,ts:10:te)

axcount = Axis(f[3,1:4],yaxisposition =:right,xticksvisible = false,yticklabelcolor = :green,ylabel="Genes in Cluster")
axn = Axis(f[3,1:4],yticklabelcolor = :blue,ylabel="Number of Clusters",xlabel="Time")

hidexdecorations!(axcount)

lines!(axcount,T,cluster_count.(d["CorrCluster"][T.+1]),color=:green)
lines!(axn,T,n_cluster.(d["CorrCluster"][T.+1]),color=:blue)

xlims!(axcount,(ts,te))
xlims!(axn,(ts,te))
ylims!(axcount,low=0)
ylims!(axn,low=0)

axn.xticks = (
    range(ts, te; length = 5),
    [string.(round.(Integer, x / 1000)) .* "K" for x in range(ts,te;length=5)]
     )


#---
# axclust = Axis(f[3,1:4])
#
# maxclust = size(d["ClusterCountHealthy"])[1]
#
# colors = vcat([:blue,:green,:red],fill(:orange,maxclust-3))
#
# for n in 1:N
#     for (i,vec) in enumerate(T[n])
#         !isempty(vec) && scatter!(
#             axclust,vec,fill(n,length(vec)),
#             color = colors[i],
#             marker = Rect, markersize = 1,
#             markerspace = :data,
#             align=(:center,:center)
#             )
#     end
# end
#
# xlims!(axclust,(ts,te))
# linkxaxes!(axprev,axclust)
# ylims!(0.5,N+0.5)

#T = get_scatter_data(d)

#

#---
# t0 = 25000
# t1 = 75000
#
# push!(colors,:gray)
#
# for (i,t) in enumerate([t0,t1,t0,t1])
#     if iseven(i)
#         label="Ill"
#         abs = 2*d["Ill"][t]
#     else
#         label="Healthy"
#         abs = 2*(d["PopSize"][t]-d["Ill"][t])
#     end
#     data = vcat(d["ClusterCount" * label][:,t],[abs-sum(d["ClusterCount" * label][:,t])])
#     ax, p = pie(
#         f[4,i],
#         data;
#         color=colors,axis=(autolimitaspect=1,),
#         axis = (; title="t=$t",xlabel=label)
#     )

    # hidexdecorations!(ax,label=false)
    # hideydecorations!(ax)
    # hidespines!(ax)
# end

# vlines!(axload,[t0,t1])

f

# save(abs_path * "/AlleleFrequenciesandProportions.pdf", f, pt_per_unit=2)

# freqs = [
#     mean(loadpos[2,i] for i in filter(!iszero,d["CorrCluster"][end][:,j]))
#     for j in 1:size(d["CorrCluster"][end])[2]
#     ]
#
# ns = [
#     length(filter(!iszero,d["CorrCluster"][end][:,j]))
#     for j in 1:size(d["CorrCluster"][end])[2]
#     ]
#
# function Prev(μ,N,fs,ns)
#     k = N-sum(ns)
#     return 1 - (1-h(N,μ)^2)^k * prod((1-fs[i]^2)^ns[i] for i in 1:length(fs))
# end
