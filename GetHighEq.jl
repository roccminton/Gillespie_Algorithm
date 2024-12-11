include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")

import .ToolkitNRM
import .PlotFromDicts

using Plots
using Measures
using CSV
using DataFrames
using Statistics

preveq(dni) = 1-exp(-dni/2)

function dosimulation(
	dnis,
	abs_path="/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/HighPrevLevel/",
	nruns=3,K=10_000,tend=100_000
	)
	#save local constants
	tspan = 10_000
	σ = 0.1

	#create new empy data frame
	mean_prev = DataFrame(dni=Float64[],N=Int64[],meanprev=Float64[])

	for dni in dnis
		Ns = round(Integer,max(11-dni*10,3))*10:10:100
		for N in Ns
			println("Currently at dni = $(round(dni,digits=2)), N = $N")
			means = Float64[]
			plots = Plots.Plot[]
			hs = DataFrame()
			for i in 1:nruns
				#Run Simulation
				h = ToolkitNRM.execute_cont(K,dni,N,tend)
				#Calculate Prev
				prev =  h.mlp["Ill"] ./ h.mlp["PopSize"]
				#find time of increase
				t_switch = findfirst(x->x>preveq(dni)+σ, prev)
				#Save Data
				hs[:,Symbol("PopSize_$i")] = h.mlp["PopSize"]
				hs[:,Symbol("Prev_$i")] = prev
				hs[:,Symbol("Load_$i")] = h.mlp["ML"] ./ h.mlp["PopSize"]
				#Create Plot
				p = PlotFromDicts.plot_MLP(h.mlp)
				#If switch occured calculate mean and add time line
				if !isnothing(t_switch)
					m = mean(view(prev,t_switch:tend))
					push!(means,m);
					vline!(p,[t_switch],label="")
					hline!(p,[m],label="")
				end
				push!(plots,p)
			end
			#calculate mean over nruns simulations
			!isempty(means) && push!(mean_prev,[dni,N,mean(means)])
			#save data and plot
			CSV.write(abs_path * "data_dni=$(round(dni,digits=2)),N=$N.csv",hs)
			savefig(plot(plots...,layout=(nruns,1),size=(600,300*nruns)),abs_path * "plot_dni=$(round(dni,digits=2)),N=$N.pdf")
		end
	end
	#save global data
	if length(dnis) == 1
		filename = "dni=$(round(dnis[1],digits=2))"
	else
		filename = "dnis=$(round(dnis[1],digits=2))-$(round(dnis[end],digits=2))"
	end
	CSV.write(
			abs_path * "Mean_Prev_High_" * filename * ".csv",
			mean_prev
	)

	return mean_prev
end
