include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")

import .ToolkitNRM
import .PlotFromDicts

using CairoMakie
using Statistics
using Distributions
using SparseArrays
using JLD

tend = 100_000
K = 10_000
N = 90
dni = 0.4

#Execute simulations
#---
#h0 = ToolkitNRM.execute_cont(K,dni,N,tend,0)
#h1 = ToolkitNRM.execute_cont(K,dni,N,tend,1)

#Save simulations
#---
abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N"
#abs_path = mkpath(abs_path)
#save(abs_path * "/r=0.jld",ToolkitNRM.convertforsaving(h0))
#save(abs_path * "/r=1.jld",ToolkitNRM.convertforsaving(h1))

#Load Data
#---
d1 = load(abs_path * "/r=1.jld")
d0 = load(abs_path * "/haploidloadhist1.jld")
#d0 = load(abs_path * "/r=0.jld")

#---
#Find extinction time
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

#mean over time
mot(data, time) = vec(mean(data[:, time], dims = 2))
#calculate load positions and average over time
mot_loadpos(data, time, N) =
    vec(mean(data[1:N, time] .+ data[N+1:2N, time], dims = 2))

#---
t_ylim = 50000
#index of extinction time for second average period
avi = 12
av = extimes[avi]:tend

#---

extimes = findextimes(d0)

prev1 = PlotFromDicts.replace_NaN(d1["Ill"] ./ d1["PopSize"])
prev0 = PlotFromDicts.replace_NaN(d0["Ill"] ./ d0["PopSize"])

load1 = PlotFromDicts.replace_NaN(d1["ML"] ./ d1["PopSize"])
load0 = PlotFromDicts.replace_NaN(d0["ML"] ./ d0["PopSize"])

#haploid load clas histograms
hist01 = mot(d0["HaploidLoadHist"], 1:extimes[1]) ./ 2K
hist02 = mot(d0["HaploidLoadHist"], av) ./ 2K
μ = sum((n - 1) * p / 2K for (n, p) in enumerate(hist02))
σ = sqrt(sum((n - μ)^2 * p / 2K for (n, p) in enumerate(hist02)))
hist1 = mot(d1["HaploidLoadHist"], 1:tend) ./ 2K

#mutation load histograms
mlhealthy1 = mot(d0["LoadHistHealthy"], 1:extimes[1]) ./ K
mlill1 = mot(d0["LoadHistIll"], 1:extimes[1]) ./ K
mlhealthy2 = mot(d0["LoadHistHealthy"], av) ./ K
mlill2 = mot(d0["LoadHistIll"], av) ./ K

#load position histogram
meanload1 = mean(d0["ML"][1:extimes[1]])
meanload2 = mean(d0["ML"][av])

loadposh1 = mot_loadpos(d0["LoadPosHealthy"], 1:extimes[1], N) ./ meanload1
loadposi1 = mot_loadpos(d0["LoadPosIll"], 1:extimes[1], N) ./ meanload1
loadposh2 = mot_loadpos(d0["LoadPosHealthy"], av, N) ./ meanload2
loadposi2 = mot_loadpos(d0["LoadPosIll"], av, N) ./ meanload2

#---
#create figure

f = Figure(resolution = (1200, 1000))

Label(
    f[0, 1:2],
    text = "No Recombination",
    font = :bold,
    tellwidth = false,
    fontsize = 30,
)
Label(
    f[0, 3],
    text = "Full Recombination",
    font = :bold,
    tellwidth = false,
    fontsize = 30,
)

#create two axis which share a x-axis
axprev0 = Axis(
    f[1, 1:2],
    yticklabelcolor = :orange,
    yaxisposition = :right,
    xaxisposition = :top,
    xticks = (
        [
            mean(extimes[1:1]),
            mean(extimes[7:7]),
            mean(extimes[9:11]),
            extimes[12:end]...,
        ],
        ["C₀-C₅", "C₆-C₇", "C₈-C₁₀", "C₁₁"],
    ),
    xticksvisible = false,
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    #title = "No Recombination"
)
axload0 = Axis(
    f[1, 1:2],
    yticklabelcolor = :red,
    ylabel = "Mutation Load",
    xlabel = "Time",
    xticks = range(0, t_ylim; length = 3),
    xtickformat = x -> string.(round.(Integer, x / 1000)) .* "K",
)

axprev1 = Axis(
    f[1, 3],
    yticklabelcolor = :orange,
    yaxisposition = :right,
    ylabel = "Prevalence", #title = "Full Recombination",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
)
axload1 = Axis(
    f[1, 3],
    yticklabelcolor = :red,
    xlabel = "Time",
    xticks = range(0, t_ylim; length = 3),
    xtickformat = x -> string.(round.(Integer, x / 1000)) .* "K",
)

hidespines!(axprev1)
hidexdecorations!(axprev1)

#set the xlims
xlims!(axload0, (0, t_ylim))
xlims!(axprev0, (0, t_ylim))
xlims!(axload1, (0, t_ylim))
xlims!(axprev1, (0, t_ylim))
ylims!(axload0, low = 0)
ylims!(axprev0, low = 0)
ylims!(axprev1, low = 0)
ylims!(axload1, low = 0)

linkyaxes!(axload0, axload1)
linkyaxes!(axprev0, axprev1)

#add extimes
vlines!(axprev0, extimes, color = :black)

#plot the data
l4 = lines!(
    axprev0,
    0:t_ylim,
    prev0[1:t_ylim+1],
    color = :orange,
    label = "Prevalence",
)
l5 = lines!(
    axload0,
    0:t_ylim,
    load0[1:t_ylim+1],
    color = :red,
    label = "Mutation Load",
)
lines!(axprev1, 0:t_ylim, prev1[1:t_ylim+1], color = :orange)
lines!(axload1, 0:t_ylim, load1[1:t_ylim+1], color = :red)

#final tweaks
hideydecorations!(axprev0, grid = true)
hideydecorations!(axload1, grid = true)

#---
#Add Haploid Load Classes Plot

Label(
    f[2, 1],
    text = L"Average over $t \in [0,\dagger_{0}]$",
    tellwidth = false,
)
Label(
    f[2, 2],
    text = L"Average over $t \in [\dagger_{%$(avi-1)},t_\mathrm{end}]$",
    tellwidth = false,
)
Label(
    f[2, 3],
    text = L"Average over $t \in [0,t_\mathrm{end}]$",
    tellwidth = false,
)

axbar01 = Axis(
    f[3, 1],#title=L"Average over $t \in [0,\dagger_{C_0}]$",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    #aspect = AxisAspect(1),alignmode = Inside()
)
axbar02 = Axis(
    f[3, 2],#title=L"Average over $t \in [\dagger_{C_{%$(avi-1)}},t_\mathrm{end}]$",
    #xlabel="Haploid Load Class",
    #aspect = AxisAspect(1)
)
axbar1 = Axis(
    f[3, 3],#title=L"Average over $t \in [0,t_\mathrm{end}]$",
    #aspect = AxisAspect(1)
)

Label(f[4, :], text = "Haploid Load Classes", font = :bold, tellwidth = false)

p(dni, N) = sqrt(1 - exp(-dni / (2N)))

#plot histograms and outlines
barplot!(axbar01, 0:N, hist01)
l1 = lines!(
    axbar01,
    0:N,
    x -> pdf(Binomial(N, p(dni, N)), x),
    color = :green,
    label = L"\mathrm{Bin}(N,h)",
)
barplot!(axbar02, 0:N, hist02)
#lines!(axbar02,0:N,x->2K*pdf(Normal(μ,σ),x),color=:red,label=L"\mathrm{Normal}(\hat{\mu},\hat{\sigma}^2)")
barplot!(axbar1, 0:N, hist1)
lines!(
    axbar1,
    0:N,
    x -> pdf(Binomial(N, p(dni, N)), x),
    color = :green,
    label = L"\mathrm{Bin}(N,h)",
)

#set lims
linkyaxes!(axbar01, axbar02, axbar1)
ylims!(axbar01, low = 0)
ylims!(axbar02, low = 0)
ylims!(axbar1, low = 0)
xlim = ceil(Integer, N / 2)
xlims!(axbar01, high = xlim)
xlims!(axbar02, high = xlim)
xlims!(axbar1, high = xlim)

#add legends
#axislegend(axbar01)
#axislegend(axbar02,position=:lt)
#axislegend(axbar1)

#final tweaks
hideydecorations!(axbar02, grid = false)
hideydecorations!(axbar1, grid = false)

#---
#Add Mutation Load Classes Plot

axml1 = Axis(
    f[5, 1],#title=L"Average over $t \in [0,\dagger_{C_0}]$",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    #aspect = AxisAspect(1)
)
axml2 = Axis(
    f[5, 2],#title=L"Average over $t \in [\dagger_{C_{%$(avi-1)}},t_\mathrm{end}]$",
    #xlabel="Mutation Load Class",
    #aspect = AxisAspect(1)
)

Label(
    f[6, 1:2],
    text = "Mutation Load Classes",
    font = :bold,
    tellwidth = false,
)

yax = 0:length(mlill1)-1

#plot histograms and outlines
l2 = barplot!(
    axml1,
    yax,
    mlill1 .+ mlhealthy1,
    color = :darkblue,
    label = "Healthy",
)
l3 = barplot!(axml1, yax, mlill1, color = :salmon, label = "Ill")
barplot!(axml2, yax, mlill2 .+ mlhealthy2, color = :darkblue, label = "Healthy")
barplot!(axml2, yax, mlill2, color = :salmon, label = "Ill")

#set lims
linkyaxes!(axml1, axml2)
ylims!(axml1, low = 0)
ylims!(axml2, low = 0)
xlim = ceil(Integer, 2N / 3)
xlims!(axml1, high = xlim)
xlims!(axml2, high = xlim)

#add legends
#axislegend(axml1)
#axislegend(axbar02,position=:lt)
#axislegend(axml2)

#final tweaks
hideydecorations!(axml2, grid = false)

#---
#Add Load Position Plot

sub_f = f[5, 3] = GridLayout()

axlp1 = Axis(
    sub_f[1, 1],
    ytickformat = x -> string.(round.(x * 100, digits = 1)) .* "%",
    title = "Mutations per Gene Positions",
)
axlp2 = Axis(
    sub_f[2, 1],
    ytickformat = x -> string.(round.(x * 100, digits = 1)) .* "%",
)

barplot!(axlp1, 1:N, loadposi1 .+ loadposh1, color = :darkblue)
barplot!(axlp1, 1:N, loadposi1, color = :salmon)
barplot!(axlp2, 1:N, loadposi2 .+ loadposh2, color = :darkblue)
barplot!(axlp2, 1:N, loadposi2, color = :salmon)

Label(
    sub_f[1, 1, Right()],
    L"before $\dagger_0$", #fontsize=18,
    rotation = pi / 2,
    padding = (3, 3, 3, 3),
)
Label(
    sub_f[2, 1, Right()],
    L"after $\dagger_{%$(avi-1)}$", #fontsize=18,
    rotation = pi / 2,
    padding = (3, 3, 3, 3),
)

ylims!(axlp1, low = 0)
linkxaxes!(axlp1, axlp2)
linkyaxes!(axlp1, axlp2)
hidexdecorations!(axlp1, grid = false)

#---
#Add Legend
Legend(
    f[7, :],
    [l4, l5, l2, l3, l1],
    ["Prevalence", "Mutation Load", "Healthy", "Ill", "Bin(N,p)"],
    tellwidth = false,
    orientation = :horizontal,
    halign = :right,
    framevisible = false,
)


#---
#Adjust Layout
rowsize!(f.layout, 1, Relative(1 / 2))
rowgap!(f.layout, 5)
colgap!(f.layout, 2, 50)

f

#save(abs_path * "/CompareRec01_Plot.pdf", f, pt_per_unit=2)
