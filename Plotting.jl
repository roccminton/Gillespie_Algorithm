"""
This module *PlotFromDicts* provides some plotting functions that take as an input
a dictionary where the keys are the trait values and the values their respective
evolutionary history in terms of subpopulation size.
"""
module PlotFromDicts

#include("DiploidModel.jl")
#import .DiploidModel

using Plots
using Measures
using Statistics
using DataFrames
using SparseArrays

export plotTSS

function plotTSS(history,time=0:length(collect(values(history))[1])-1)
	p = plot(legend = false)
	min_x,max_x = min(keys(history)...),max(keys(history)...)

	color_gradient = cgrad(:thermal)
	ncolors = length(color_gradient)-1

	c(x,min,max,N) = min == max ? 1 : floor(Integer,((x - min) / (max - min)) * N)+1

	for (x, his_x) in history
		plot!(p,time,his_x,color=color_gradient.colors.colors[c(x,min_x,max_x,ncolors)])
	end

	return p
end

function plotMLPoverX(ML,P,X)
	p = plot(size = (600,300),
			rightmargin = 20mm,
			leftmargin = 15mm,
			topmargin = 5mm,
			bottommargin = 5mm,
			title = "",
			titlefontsize = 8,
			titleposition=:left
			)
	plotMLPoverX!(p,ML,P,X)
	return p
end

function plotMLPoverX!(p,ML,P,X)
	plot!(p,X,P,
		label = "", grid = false,
		color = :orange,
		ylabel = "Prevalence",
		ylim = (0,maximum(P)*1.01),
		framestyle = :zerolines,
	   )
	#plot mutation load
	maxload = ceil(Integer,maximum(ML))
	plot!(twinx(),X,ML,
		label="",grid = false,
		color=:red,
		#ymirror = false,
		ylabel = "Mutation Load",
		ylim = (0,ceil(Integer,maximum(ML))),
		xticks = false,
		framestyle = :zerolines,
		showaxis = false
	   )
	 vline!(p,[X[end]],linewidth=1,color=:black,label="")
end

function plotmutationloadandprevalence(popsize,prevalence,mutationload,time=0:(length(popsize)-1))
	#setting values
	ticksfontsize = 5

	#setup plot
	p = plot(size = (600,300),
			rightmargin = 20mm,#20mm,
			leftmargin = 15mm,#15mm,
			#topmargin = 5mm,
			bottommargin = 5mm,
			title = "",
			titlefontsize = 5,
			titleposition=:left,
			xtickfontsize=ticksfontsize,
			)
	#plot population size in background
	popticks = floor.(Int,eventicks(time[end],4))
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
	prevmax = ceil(Integer,maximum(prevalence)*100)/100
	prevticks = eventicks(prevmax,5)
	plot!(twinx(),prevalence,
		label = "", grid = false,
		color = :orange,
		ylabel = "Prevalence",
		#ylim = (0,prevmax*1.01),
		yticks = (
				prevticks,
				string.(round.(100 .* prevticks ,digits=2)).*"%"
			),
		ytickfontsize=ticksfontsize,
		#tickfontcolor = :orange,
		framestyle = :zerolines,
		xticks = false,
		showaxis = false
 	   )
	#plot mutation load
	maxload = ceil(Integer,maximum(mutationload))
	loadticks = eventicks(maxload,5)
	plot!(twinx(),mutationload,
		label="",grid = false,
		color=:red,
		ymirror = false,
		ylabel = "Mutation Load",
		#ylim = (0,maxload),
		yticks = loadticks,
		#tickfontcolor = :red,
		ytickfontsize=ticksfontsize,
		xticks = false,
		framestyle = :zerolines,
		showaxis = false
 	   )
	vline!([time[end]],linewidth=1,color=:black,label="")
end

function plotmutationloadandprevalence(history,time=0:length(collect(values(history))[1])-1)
	plotmutationloadandprevalence(
		DiploidModel.populationsize(history),
		DiploidModel.rescale(DiploidModel.abs_ill_individual(history),popsize),
		DiploidModel.rescale(DiploidModel.abs_mutationload(history),popsize),
		time)
end

function mutationdistributiongif(history)
	#clean data
	total_mutations_per_time_rescaled = DiploidModel.rel_mutation_distribution(history)
	populationsize = DiploidModel.populationsize(history)
	#set maximum
	m = ceil(maximum(maximum.(total_mutations_per_time_rescaled)),digits=2)
	nLoci = length(first(first(history)).genes)

	#resized time
	time = range(1,nLoci; length=length(populationsize))

	anim = @animate for (i,mutations) in enumerate(total_mutations_per_time_rescaled)
		#plot population size in background
		plot(
			time,
			populationsize,
			legend=false,
			color = :gray,
			xticks = false,
			yticks = false,
			leftmargin = 10mm,
			bottommargin = 10mm,
			)
		#add scatter to current population size
		scatter!([time[i]],[populationsize[i]],color=:red)
		#add scatter plot of gene distribution
		bar!(
			twinx(),mutations,
			ylim=(0,m),ymirror = false,
			label="Load per Locus",
			xlabel = "Loci"
			)
		#add mean mutation load per gene
		hline!([mean(mutations)],color=:red,label="mean")
	end
	return anim
end

function eventicks(stop,nticks,start=0)
	l = stop-start
	step = round(l/nticks,digits=-(floor(Int,log10(l)-1)))
	return range(start,stop;step = step)
end

"""
Generates the Mutation Load and Prevalence Plot.
"""
plot_MLP(history) = plotmutationloadandprevalence(
                history["PopSize"],
                replace_NaN(history["Ill"] ./ history["PopSize"]),
                replace_NaN(history["ML"] ./ history["PopSize"]))
plot_MLP(df::DataFrame) = plotmutationloadandprevalence(
                df.PopSize,
                replace_NaN(df.Ill ./ df.PopSize),
                replace_NaN(df.Mutation ./ df.PopSize))
function plot_MLP(history,abs_path)
	p = plot_MLP(history)
	savefig(p,abs_path)
end

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)

function plot_MLHist(hhist,t,popsize,ylim,xlim,key="both")
	if key == "both"
		h = hhist["Ill"][t] .+ hhist["Healthy"][t]
		color = :blue
	else
		h = hhist[key][t]
		color = key=="Ill" ? :orange : :red
	end
	class, abs_freq = findnz(h)
	#Indexing statrs with 1 but mutation load starts with 0
	class .= class .- 1
	rel_freq = abs_freq ./ popsize
	bar(
		class, rel_freq,
		xlabel="Classes",
		ylabel="Frequency",
		legend=false,
		color = color,
		xlim = (-1,xlim),
		ylim = (0,ylim),
	)
end

function gif_MLP_LoadHist(history,abs_path,tend;everyn=100,maxfreq=0.25)
	maxmut = maxmutatioins(history.loadhist)
	anim = @animate for t in 1:tend
	    p1 = plot_MLP(history.mlp)
	    vline!([t])
	    p2 = plot_MLHist(history.loadhist,t,history.mlp["PopSize"][t],maxfreq,maxmut)
	    vline!([history.mlp["ML"][t]/history.mlp["PopSize"][t]])
	    plot(p1,p2,layout=(2,1),legend=false)
	end every everyn

	gif(anim, abs_path)
end

function plot_LoadPos(loadpos,t,popsize,ylim,xlim;key="both",orientation=:v)
	if key == "both"
		h = loadpos["Ill"][t] .+ loadpos["Healthy"][t]
		color = :blue
	else
		h = loadpos[key][t]
		color = key=="Ill" ? :orange : :red
	end
	if orientation == :h
		xlim,ylim = ylim,xlim
		xlabel, ylabel = "Frequency", "Loci"
	else
		xlabel, ylabel = "Loci", "Frequency"
	end
	p = plot(
		ylim=ylim,xlim=xlim,
		legend=false,
		xlabel=xlabel,ylabel=ylabel,
		)
	plot_LoadPos!(p,h,popsize,color,orientation)
	return p
end

function plot_LoadPos!(p,h::SparseVector,ps,c,i,orientation)
	class, abs_freq = findnz(h)
	rel_freq = i .* (abs_freq ./ ps)
	bar!(p,
		class, rel_freq,
		label="",
		color = c,
		orientation = orientation
	)
end

function plot_LoadPos!(p,h::Vector{SparseVector{S,T}},ps,c,orientation) where {S<:Number,T<:Number}
	plot_LoadPos!(p,h[1],ps,c,1,orientation)
	plot_LoadPos!(p,h[2],ps,c,-1,orientation)
end

function gif_MLP_LoadPos(history,abs_path,tend;everyn=100,maxfreq=0.5)
	anim = @animate for t in 1:tend
	    p1 = plot_MLP(history.mlp)
	    vline!([t])
	    p2 = plot_LoadPos(
			history.loadpos,t,
			history.mlp["PopSize"][t],
			(-maxfreq,maxfreq),
			(0,history.par.Nloci+1)
			)
	    plot(p1,p2,layout=(2,1),legend=false)
	end every everyn

	gif(anim, abs_path)
end

function gif_MLP_LoadPos_LoadHis(history,abs_path,tend;everyn=100,maxfreqhis=0.25,maxfreqpos=0.5)

	maxmut = maxmutatioins(history.loadhist)

	anim = @animate for t in 1:tend
		pos = plot_LoadPos(
		    history.loadpos,
		    t,
		    history.mlp["PopSize"][t],
		    (-maxfreqpos, maxfreqpos),
		    (0, history.par.Nloci+1),
		    orientation = :h,
		)
		his = plot_MLHist(
		    history.loadhist,
		    t,
		    history.mlp["PopSize"][t],
		    maxfreqhis,maxmut,
		    )
		mlp = plot_MLP(history.mlp)
		vline!([t],label="")

		l = @layout [[a ; b] c{0.3w} ]

		plot(
			mlp,his,pos,
			layout=l,
			size=(800,400)
			)
	end every everyn

	gif(anim, abs_path)
end

function maxmutatioins(hhist,key="both")
    if key=="both"
        return max(findnz(hhist["Ill"][end])[1][end],findnz(hhist["Healthy"][end])[1][end])
    else
        return findnz(hhist[key][end])[1][end]
    end
end

end #end of Module PlotFromDicts
