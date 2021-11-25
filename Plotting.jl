"""
This module *PlotFromDicts* provides some plotting functions that take as an input
a dictionary where the keys are the trait values and the values their respective
evolutionary history in terms of subpopulation size.
"""
module PlotFromDicts

using Plots

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

end #end of Module PlotFromDicts
