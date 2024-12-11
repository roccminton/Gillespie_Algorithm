include("NonRandomMating/Toolkit.jl")

import .ToolkitNRM

using CairoMakie
using Statistics
using Distributions
using SparseArrays
using JLD
using LinearAlgebra

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
#d = load(abs_path * "/alldata1.jld")

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

cluster_count(s) = iszero(length(s)) ? 0 : sum(count(!iszero, s, dims = 1))

#---
extimes = findextimes(d)
#mutinctimes = [69000,72000,90000]

times = [1:extimes[1], extimes[1]:extimes[9], extimes[end]:tend]
#times = [extimes[6]:extimes[7],extimes[9]:20000]
# t0,t1,t2 = 38_500, 43_500,65_000
# times = [t0:t1,t1:t2]

#start and end of plot window
ts = 0
te = tend

#if compareplot should be printed
compl = false

#---

prev = replace_NaN(d["Ill"] ./ d["PopSize"])
ml = replace_NaN(d["ML"] ./ d["PopSize"])

cormatrixs = [CovToCor(mean(d["CovMatrixList"][int])) for int in times]
compl && push!(cormatrixs, cormatrixs[2] .- cormatrixs[1])

loadpos = [mean(d["LoadPos2"][:, :, int], dims = 3) ./ K for int in times]
compl && push!(loadpos, loadpos[2] .- loadpos[1])

# loadpos_h = [
#     mot_loadpos(d["LoadPosHealthy"], times[i]:times[i+1], N) ./ K for i in 1:length(times)-1
# ]
# loadpos_i = [
#     mot_loadpos(d["LoadPosIll"], times[i]:times[i+1], N) ./ K for i in 1:length(times)-1
# ]


#---
#create figure

f = Figure()

# a = f[1,1] = GridLayout()
# c = f[2:3, 1] = GridLayout()

#---

#create two axis which share a x-axis
# axprev = Axis(
#     f[1, 1:4],
#     yticklabelcolor = :orange,
#     yaxisposition = :right,
#     xaxisposition = :top,
#     xticksvisible = false,
#     ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
#     ylabel = "Prevalence",
# )
#
# if compl
#     axprev.xticks = ([t0, t1, t2], [L"t_0", L"t_1", L"t_2"])
# else
#     axprev.xticks = (extimes, [L"\dagger" for i = 0:length(extimes)-1])
# end
#
# axload = Axis(
#     f[1, 1:4],
#     yticklabelcolor = :red,
#     ylabel = "Mutation Load",
#     xlabel = "Time",
# )
#
# if compl
#     axload.xticks = (
#         vcat(range(ts, te; length = 5), extimes),
#         vcat(
#             [
#                 string.(round.(Integer, x / 1000)) .* "K" for
#                 x in range(ts, te; length = 5)
#             ],
#             [L"\dagger" for i = 1:length(extimes)],
#         ),
#     )
# else
#     axload.xticks = (
#         range(ts, te; length = 5),
#         [
#             string.(round.(Integer, x / 1000)) .* "K" for
#             x in range(ts, te; length = 5)
#         ],
#     )
# end
#
# # hidespines!(axprev)
# # hidexdecorations!(axprev)
#
# #set the xlims
# linkxaxes!(axload, axprev)
# xlims!(axload, (ts, te))
# xlims!(axprev, (ts, te))
# !compl && ylims!(axload, low = 0)
# !compl && ylims!(axprev, low = 0)
#
# #vlines!(axload,mutinctimes,color=:darkred,linewidth=1)
#
# #plot the data
# prevline = lines!(
#     axprev,
#     ts:te,
#     prev[ts+1:te+1],
#     color = :orange,
#     label = "Prevalence",
# )
# loadline =
#     lines!(axload, ts:te, ml[ts+1:te+1], color = :red, label = "Mutation Load")


##
rows = 1
cols = 3

if compl
    titles = [
        L"Average over $t \in [t_0,t_1]$",
        L"$\leftarrow$ Difference $\rightarrow$",
        L"Average over $t \in [t_1,t_2]$",
    ]
    order = [1, 3, 2]
else
    titles = [
        L"Average over $t \in [0,\dagger_0]",
        L"Average over $t \in [\dagger_0,\dagger_{11}]",
        L"Average over $t \in [\dagger_{11},t_{\mathrm{end}}]",
    ]
    order = [1, 2, 3]
end

for (i, o) in enumerate(order)
    col = (rem(i, cols) == 0 ? cols : rem(i, cols))
    #ax = f[div(i + cols - 1, cols)+1, col] = GridLayout()
    ax = f[1, col] = GridLayout()
    colsize!(f.layout, col, Aspect(1, 1.0))
    ax1 = Axis(
        ax[1:2, 1],
        title = titles[i],
        ylabel = isone(i) ? "Correlation Matrix" : "",
    )
    heatmap!(ax1, cormatrixs[o], clim = compl ? (-1, 1) : (0, 1))
    ax2 = Axis(
        ax[3, 1],
        #xlabel = i == 2 ? "Muation Distribution per Locus" : "",
        ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    )
    hidespines!(ax2)
    xlims!(ax2, (0.5, N + 0.5))
    barplot!(
        ax2,
        loadpos[o][1, :] .+ loadpos[o][2, :] .+ loadpos[o][3, :],
        color = :gray,
        label = "0/0",
    )
    barplot!(
        ax2,
        loadpos[o][2, :] .+ loadpos[o][3, :],
        color = :darkblue,
        label = "0/1",
    )
    barplot!(ax2, loadpos[o][3, :], color = :salmon, label = "1/1")
    rowgap!(ax, 1)
    hidexdecorations!(ax1)
    hideydecorations!(ax1, label = false)
    isone(i) ? hideydecorations!(ax2, ticklabels = false, ticks = false) :
    hideydecorations!(ax2)
    hidexdecorations!(ax2, label = false)
end

Colorbar(f[1, 4], limits = compl ? (-1, 1) : (0, 1))

#---,xticksvisible = false

T = round.(Integer,ts:10:te)

axcount = Axis(
    f[3,1:4],
    xaxisposition = :top,xticksvisible = false,
    yaxisposition =:right,yticklabelcolor = :green,ylabel="Genes in Cluster"
    )
axn = Axis(f[3,1:4],yticklabelcolor = :blue,ylabel="Number of Clusters",xlabel="Time")

hidexdecorations!(axcount,ticks=false,ticklabels=false,grid=false)

#lines!(axcount,T,cluster_count.(d["CorrCluster"][T.+1]),color=:green)
#lines!(axn,T,n_cluster.(d["CorrCluster"][T.+1]),color=:blue)

xlims!(axcount,(ts,te))
xlims!(axn,(ts,te))
ylims!(axcount,low=0)
ylims!(axn,low=0)

axn.xticks = (
    range(ts, te; length = 3),
    [string.(round.(Integer, x / 1000)) .* "K" for x in range(ts,te;length=3)]
     )

axcount.xticks = (extimes, [L"\dagger" for i = 0:length(extimes)-1])

#---


m = :rect
Label(f[2, 1:2], text = "Mutation Distribution per Locus", font = :bold, tellwidth = false)
Legend(
    f[2, 3],
    [
        MarkerElement(color = :gray, marker = m),
        MarkerElement(color = :darkblue, marker = m),
        MarkerElement(color = :salmon, marker = m),
    ],
    ["0/0", "0/1", "1/1"],
    tellwidth = false,
    orientation = :horizontal,
    halign = :right,
    framevisible = false,
)


colgap!(f.layout, 10)
rowgap!(f.layout, 10)

resize_to_layout!(f)

f

#save(abs_path * "/ClusterandMatrix.pdf", f, pt_per_unit=2)
