"""
This module *PlotFromDicts* provides some plotting functions that take as an input
a dictionary where the keys are the trait values and the values their respective
evolutionary history in terms of subpopulation size.
"""
module PlotFromDicts

include("DiploidModel.jl")
import .DiploidModel

using Plots
using Measures
using Statistics

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

end #end of Module PlotFromDicts
