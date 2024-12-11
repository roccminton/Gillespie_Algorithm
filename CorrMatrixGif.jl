include("NonRandomMating/Toolkit.jl")
include("ToolBoxPlotting.jl")


import .ToolkitNRM

using CairoMakie
using Statistics
using Distributions
using SparseArrays
using JLD
using LinearAlgebra

tend = 100_000
K = 10_000
N = 400
dni = 0.07
i = 3

#Execute simulations
#---
#h = ToolkitNRM.execute_cont(K,dni,N,tend,0)

#Save simulations
#---
abs_path = "/media/larocca/PortableSSD/Data/NoRecombination/CorrMatrix/N=$N,dni=$dni"
#abs_path = mkpath(abs_path)
filename = "data_$i"
#save(abs_path * "/" * filename * ".jld",ToolkitNRM.convertforsaving(h))

#Load Data (takes up to 45 min)
#---

#d = load(abs_path * "/" * filename * ".jld")
#M = load(abs_path * "/" * filename * "CovMatrixs.jld")


#---
extimes = findextimes(d)

#start and end of plot window
ts = 14_000
te = 25_000

corrkeys = sort(parse.(Int,collect(keys(M))))

#---

prev = replace_NaN(d["Ill"] ./ d["PopSize"])
ml = replace_NaN(d["ML"] ./ d["PopSize"])

#cormatrixs = CovToCor.(d["CovMatrixList"])

t = Observable(1)

#---
#create figure
f = Figure(resolution=(600,800))
#---
#create two axis which share a x-axis
axprev = Axis(
    f[1, 1:2],
    yticklabelcolor = :orange,
    yaxisposition = :right,
    xaxisposition = :top,
    xticksvisible = false,
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    ylabel = "Prevalence",
)

axprev.xticks = (extimes, [L"\dagger" for i = 0:length(extimes)-1])

axload = Axis(
    f[1, 1:2],
    yticklabelcolor = :red,
    ylabel = "Mutation Load",
    xlabel = "Time",
)

axload.xticks = (
        range(ts, te; length = 5),
        [
            string.(round.(Integer, x / 1000)) .* "K" for
            x in range(ts, te; length = 5)
        ],
    )
# hidespines!(axprev)
# hidexdecorations!(axprev)

#set the xlims
linkxaxes!(axload, axprev)
xlims!(axload, (ts, te))
xlims!(axprev, (ts, te))

vlines!(axload,@lift([$t]),color=:black,linewidth=1)

#plot the data
prevline = lines!(axprev, ts:te, prev[ts+1:te+1],color = :orange,label = "Prevalence")
loadline = lines!(axload, ts:te, ml[ts+1:te+1], color = :red, label = "Mutation Load")

##

axcor = Axis(f[2,1:2])
heatmap!(axcor,@lift(cormatrixs[$t]),clim=(0,1))

hidedecorations!(axcor)

Colorbar(f[2, 2], limits = (0, 1))

#---

T = ts:te

axcount = Axis(
    f[3,1:2],
    xaxisposition = :top,xticksvisible = false,
    yaxisposition =:right,yticklabelcolor = :green,ylabel="Genes in Cluster"
    )
axn = Axis(f[3,1:2],yticklabelcolor = :blue,ylabel="Number of Clusters",xlabel="Time")

hidexdecorations!(axcount,ticks=false,ticklabels=false,grid=false)

lines!(axcount,T,cluster_count.(d["CorrCluster"][T.+1]),color=:green)
lines!(axn,T,n_cluster.(d["CorrCluster"][T.+1]),color=:blue)

xlims!(axcount,(ts,te))
xlims!(axn,(ts,te))
ylims!(axcount,low=0)
ylims!(axn,low=0)

axn.xticks = (
    range(ts, te; length = 3),
    [string.(round.(Integer, x / 1000)) .* "K" for x in range(ts,te;length=3)]
     )

axcount.xticks = (extimes, [L"\dagger" for i = 0:length(extimes)-1])

vlines!(axcount,@lift([$t]),color=:black,linewidth=1)

resize_to_layout!(f)

#---

#timestamps = round.(Integer,range(ts,te,length=500))
timestamps = ts:10:te

# record(f,abs_path * "/CovMatrixZoomSlow.gif", timestamps) do x
#     t[]=x
# end

f

#save(abs_path * "/ClusterandMatrix.pdf", f, pt_per_unit=2)
