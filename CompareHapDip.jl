include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")
include("HaighSimulations.jl")


import .ToolkitNRM
import .PlotFromDicts

using Plots
using Measures

tend = 50_000
K = 10000
N = 75
dni = 0.6

h0 = ToolkitNRM.execute_wf(K,dni,N,tend,0)
#---
s = 0.14
λ = 0.66
slots = 50

h1 = runHaigh(vcat([0,K],zeros(Integer,slots-1)),s,λ,tend)
#---

function plotmutationloadandprevalence(
	popsize,prevalence,mutationload,legend,label1,label2,
	prevmax, maxload,rm,lm,title,
	time=0:(length(popsize)-1)
	)
	#setting values
	ticksfontsize = 7#5

	#setup plot
	p = plot(size = (600,300),
			rightmargin = rm,#20mm,
			leftmargin = lm,#15mm,
			#topmargin = 5mm,
			bottommargin = 5mm,
			title = title,
			titlefontsize = 15,
			titleposition=:center,
			xtickfontsize=ticksfontsize,
			legend = :topleft,
			)
	#plot population size in background
	popticks = floor.(Int,PlotFromDicts.eventicks(time[end],2))
	plot!(p,time,popsize,
		grid = false, label = "",
		color = Gray(0.2), alpha = 0.7,
		yticks = false,
		ylims = (0,maximum(popsize)*1.01),
		xlabel = "Time",
		xticks = (popticks,string.(popticks)),
		framestyle = :zerolines,
    	)
	#plot Prevalence
	prevticks = PlotFromDicts.eventicks(prevmax,5)
	plot!(twinx(),prevalence,
		label = legend ? "Prevalence" : "",
		grid = false,
		color = :orange,
		ylabel = label2 ? "Prevalence" : "",
		ylim = (0,prevmax*1.01),
		#ylim = (0,0.60),
		yticks = label2 ? (
				prevticks,
				string.(round.(100 .* prevticks ,digits=2)).*"%"
			) : false,
		ytickfontsize=ticksfontsize,
		#tickfontcolor = :orange,
		framestyle = :zerolines,
		xticks = false,
		showaxis = false,
 	   )
	#plot mutation load
	loadticks = PlotFromDicts.eventicks(maxload,5)
	plot!(twinx(),mutationload,
		label=legend ? "Mutation Load" : "",
		grid = false,
		color=:red,
		ymirror = false,
		ylabel = label1 ? "Mutation Load" : "",
		ylim = (0,maxload),
		yticks = label1 ? loadticks : false,
		#tickfontcolor = :red,
		ytickfontsize=ticksfontsize,
		xticks = false,
		framestyle = :zerolines,
		showaxis = false
 	   )
	vline!([time[end]],linewidth=1,color=:black,label="")
end

#----

prev1 = 1 .- [mean_fitness(h1[i,:],s) for i in 1:tend]
prev0 = PlotFromDicts.replace_NaN(h0["Ill"] ./ h0["PopSize"])

load1 = [load(h1[i,:]) for i in 1:tend] ./ K
load0 = PlotFromDicts.replace_NaN(h0["ML"] ./ h0["PopSize"])

prevmax = max(ceil(Integer,maximum(prev1)*100)/100,ceil(Integer,maximum(prev0)*100)/100)
maxload = max(ceil(Integer,maximum(load1)),ceil(Integer,maximum(load0)))

p1 = plotmutationloadandprevalence(
	fill(K,tend),
	prev1,
	load1,
	false,true,false,
	prevmax, maxload,
	0mm, 15mm,
	"haploid"
	)

p0 = plotmutationloadandprevalence(
	h0["PopSize"],
	prev0,
	load0,
	false,false,true,
	prevmax, maxload,
	20mm, 0mm,
	"diploid"
	)

plot(p1,p0,layout=(1,2),size=(800,450))
savefig("/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/WF_Compare_HapDip_N=$N,dni=$dni.pdf")
