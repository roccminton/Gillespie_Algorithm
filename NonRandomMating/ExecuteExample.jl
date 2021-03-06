include("Toolkit.jl")
include("../Plotting.jl")

import .ToolkitNRM
import .PlotFromDicts

using Plots

#(constant) Population size
K = 10_000
#Number of loci/genes
N = 60
#Number of Generations
tend = 5_000
#de novo mutation rate
dni = 0.8


#---
#run one simulation
run_simulation = true

#run several simulations an take the average
run_simulations = false
#number of runs to average over
n = 50

#save result in file
save_data = false
#path to save data at
abs_path = "/home/larocca/github/Gillespie_Algorithm/NonRandomMating/Output/ConstPopSize/"

#plot result
plot_data = true

#---

#run simulations with mean
if run_simulations
    wf_history = ToolkitNRM.execute_and_mean(K,dni,N,tend,n,ToolkitNRM.execute_wf)
    cont_history = ToolkitNRM.execute_and_mean(K,dni,N,tend,n,ToolkitNRM.execute_cont)
end

#run only one simulation each
if run_simulation
    wf_history = ToolkitNRM.execute_wf(K,dni,N,tend)
    cont_history = ToolkitNRM.execute_cont(K,dni,N,tend)
end

#save data
if save_data
    ToolkitNRM.save_data(abs_path,"Cont_K=$K,dni=$dni,Nloci=$N",cont_history,tend)
    ToolkitNRM.save_data(abs_path,"WF_K=$K,dni=$dni,Nloci=$N",wf_history,tend)
end

#plot data
if plot_data
    plot(ToolkitNRM.plot_MLP(wf_history),ToolkitNRM.plot_MLP(cont_history),layout=(2,1))
end
