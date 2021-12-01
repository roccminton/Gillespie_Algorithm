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

function plotmutationloadandprevalence(history,time=0:length(collect(values(history))[1])-1)
	popsize = DiploidModel.populationsize(history)
	#setup plot
	p = plot(size = (600,300),
			rightmargin = 20mm,
			leftmargin = 15mm,
			topmargin = 5mm,
			bottommargin = 5mm,
			title = "",
			titlefontsize = 8,
			titleposition=:left
			)
	#plot population size in background
	plot!(p,time,popsize,
		grid = false, label = "",
		color = Gray(0.2), alpha = 0.7,
		yticks = false,
		ylims = (0,maximum(popsize)*1.01),
		xlabel = "Time",
		xticks = floor.(Int,range(0,time[end];length=5)),
		framestyle = :zerolines,
    	)
	#plot Prevalence
	prevalence = DiploidModel.ill_individual(history) ./ popsize
	prevmax = ceil(Integer,maximum(prevalence)*100)/100
	plot!(twinx(),prevalence,
		label = "", grid = false,
		color = :orange,
		ylabel = "Prevalence",
		ylim = (0,maximum(prevalence)*1.01),
		yticks = (
				range(0,prevmax;length=5),
				string.(round.(100 .*range(0,prevmax;length=5),digits=2)).*"%"
			),
		#tickfontcolor = :orange,
		framestyle = :zerolines,
		xticks = false,
		showaxis = false
 	   )
	#plot mutation load
	mutationload = DiploidModel.mutationload(history) ./ popsize
	maxload = ceil(Integer,maximum(mutationload))
	plot!(twinx(),mutationload,
		label="",grid = false,
		color=:red,
		ymirror = false,
		ylabel = "Mutation Load",
		ylim = (0,maxload),
		yticks = range(0,maxload;length=5),
		#tickfontcolor = :red,
		xticks = false,
		framestyle = :zerolines,
		showaxis = false
 	   )
	vline!([time[end]],linewidth=1,color=:black,label="")
end

end #end of Module PlotFromDicts
