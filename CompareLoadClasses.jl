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

d = load(abs_path * "/corrcluster1.jld")

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
t_ylim = tend
#index of extinction time for second average period
extimes = findextimes(d)
avi = length(extimes)
av = extimes[end]:tend

#---



prev = PlotFromDicts.replace_NaN(d["Ill"] ./ d["PopSize"])
ml = PlotFromDicts.replace_NaN(d["ML"] ./ d["PopSize"])

#haploid  clas histograms
hist1 = mot(d["HaploidLoadHist"], 1:extimes[1]) ./ 2K
hist2 = mot(d["HaploidLoadHist"], av) ./ 2K
μ = sum((n - 1) * p / 2K for (n, p) in enumerate(hist2))
σ = sqrt(sum((n - μ)^2 * p / 2K for (n, p) in enumerate(hist2)))


#mutation load histograms
mlhealthy1 = mot(d["LoadHistHealthy"], 1:extimes[1]) ./ K
mlill1 = mot(d["LoadHistIll"], 1:extimes[1]) ./ K
mlhealthy2 = mot(d["LoadHistHealthy"], av) ./ K
mlill2 = mot(d["LoadHistIll"], av) ./ K

#load position histogram
meanload1 = mean(d["ML"][1:extimes[1]])
meanload2 = mean(d["ML"][av])

loadposh1 = mot_loadpos(d["LoadPosHealthy"], 1:extimes[1], N) ./ meanload1
loadposi1 = mot_loadpos(d["LoadPosIll"], 1:extimes[1], N) ./ meanload1
loadposh2 = mot_loadpos(d["LoadPosHealthy"], av, N) ./ meanload2
loadposi2 = mot_loadpos(d["LoadPosIll"], av, N) ./ meanload2

#---
#create figure

f = Figure(resolution = (1200, 1000))

#create two axis which share a x-axis
axprev = Axis(
    f[1, 1:2],
    yticklabelcolor = :orange,
    yaxisposition = :right,
    xaxisposition = :top,
    xticks = (extimes,fill(L"\dagger",length(extimes))),
    xticksvisible = false,
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    ylabel = "Prevalence"
)
axload = Axis(
    f[1, 1:2],
    yticklabelcolor = :red,
    ylabel = "Mutation Load",
    xlabel = "Time",
    xticks = range(0, t_ylim; length = 3),
    xtickformat = x -> string.(round.(Integer, x / 1000)) .* "K",
)

#set the xlims
xlims!(axload, (0, t_ylim))
xlims!(axprev, (0, t_ylim))
ylims!(axload, low = 0)
ylims!(axprev, low = 0)

#plot the data
l4 = lines!(
    axprev,
    0:t_ylim,
    prev[1:t_ylim+1],
    color = :orange,
    label = "Prevalence",
)
l5 = lines!(
    axload,
    0:t_ylim,
    ml[1:t_ylim+1],
    color = :red,
    label = "Mutation Load",
)

eq_load = 2*N*sqrt(1-exp(-dni/(2N)))
eq_prev = 1-exp(-dni/2)

hlines!(axload,[eq_load],color=:darkred)
hlines!(axprev,[1-exp(-dni/2)],color=:darkorange)

text!(
    axload,extimes[1],eq_load,
    text=L"\hat{L} = 2N\sqrt{1-e^{-\frac{\mu}{N}}}",
    color=:darkred,offset=(20,-50),
    fontsize = 30, align = (:left,:center)
    )

text!(
    axprev,extimes[1],eq_prev,
    text=L"\hat{P} = 1-e^{-\mu}",
    color=:darkorange,offset=(20,-30),
    fontsize = 30, align = (:left,:center)
    )

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

axbar1 = Axis(
    f[3, 1],#title=L"Average over $t \in [0,\dagger_{C_0}]$",
    ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
    #aspect = AxisAspect(1),alignmode = Inside()
)
axbar2 = Axis(
    f[3, 2],#title=L"Average over $t \in [\dagger_{C_{%$(avi-1)}},t_\mathrm{end}]$",
    #xlabel="Haploid Load Class",
    #aspect = AxisAspect(1)
)

Label(f[4, :], text = "Haploid Load Classes", font = :bold, tellwidth = false)

p(dni, N) = sqrt(1 - exp(-dni / (2N)))

#plot histograms and outlines
barplot!(axbar1, 0:N, hist1)
l1 = lines!(
    axbar1,
    0:N,
    x -> pdf(Binomial(N, p(dni, N)), x),
    color = :green,
    label = L"\mathrm{Bin}(N,h)",
)

text!(axbar1,N/4,0.1,text=L"p=\sqrt{1-e^{-\frac{\mu}{N}}}",color=:green,fontsize=30,align = (:left,:center))
barplot!(axbar2, 0:N, hist2)
# lines!(axbar2,0:N,x->2K*pdf(Normal(μ,σ),x),color=:red,label=L"\mathrm{Normal}(\hat{\mu},\hat{\sigma}^2)")

#set lims
linkyaxes!(axbar1, axbar2)
ylims!(axbar1, low = 0)
ylims!(axbar2, low = 0)
linkxaxes!(axbar1,axbar2)
xlim = ceil(Integer,N/2)
xlims!(axbar1, (-0.5,xlim+0.5))
xlims!(axbar2, (-0.5,xlim+0.5))

hideydecorations!(axbar2, grid = false)

#add legends
# axislegend(axbar2,position=:lt)
# axislegend(axbar1)


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
    f[6, :],
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
linkxaxes!(axml1,axml2)
xlims!(axml1, low = -0.5)
xlims!(axml2, low = -0.5)

#add legends
# axislegend(axml1)
# axislegend(axbar2,position=:lt)
# axislegend(axml2)

#final tweaks
hideydecorations!(axml2, grid = false)


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
# colgap!(f.layout, 2, 50)

f

#save(abs_path * "/LoadClassComp.pdf", f, pt_per_unit=2)
