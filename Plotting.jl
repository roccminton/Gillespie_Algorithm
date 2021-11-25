
module PlotFromDicts

using Plots

export plotTSS

function plotTSS(history,time=0:length(collect(values(history))[1])-1)
	p = plot(legend = false)
	min_x,max_x = min(keys(history)...),max(keys(history)...)

	color_gradient = cgrad(:thermal)
	ncolors = length(color_gradient)-1

	c(x) = floor(Integer,((x - min_x) / (max_x - min_x)) * ncolors)+1

	for (x, his_x) in history
		plot!(p,time,his_x,color=color_gradient.colors.colors[c(x)])
	end

	return p
end

end #end of Module PlotFromDicts
