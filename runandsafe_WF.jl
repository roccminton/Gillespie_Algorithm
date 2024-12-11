include("WrightFisher.jl")
import .WrightFisher

#choose birht death functions package
include("/home/larocca/github/Gillespie_Algorithm/WrightFisher_Diploid_NoRec_Independent.jl")
#choose stats functions package
include("/home/larocca/github/Gillespie_Algorithm/WF_Stats/HLC.jl")
#plotting functions
include("/home/larocca/github/Gillespie_Algorithm/ToolBoxPlotting.jl")

#---
using JLD
using CairoMakie

#---

model_parameter(dni,N,K) = (
        μ=dni, Nloci=N, K=K,
        popsize=const_popsize,
        )

function convertforsaving(h,N,dni,K)
    return merge(
            h,
            Dict("μ" => dni, "K" => K, "Nloci" => N)
        )
end

#---
function dosimulation(
        N,dni,filename;
        K=10_000,tend=100_000,nruns=3,
        abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/"
        )

    abs_path = mkpath(abs_path * "K=$K,dni=$dni,N=$N")
    filename = "/" * filename

    println("Currently at dni = $(round(dni,digits=2)), N = $N")
    f = Figure(resolution=(600,300*nruns))
    for i in 1:nruns
        #Run Simulation
        n0 = healthy_pop(K)
        h = WrightFisher.run_wrightfisher(
            1:tend+1,
            n0,
            model_parameter(dni,N,K)
        )
        #Add cumultative statistics
        add_stats!(h,n0)
        #Save Data
        save(abs_path * filename *"$i.jld",convertforsaving(h,N,dni,K))
        #Create Plot
        add_mlp_plot!(f,h,tend,i=i)
    end
    #save data and plot
    save(abs_path * filename *"_overview.pdf", f)
end

function dolongsimulation(
        N,dni,filename;
        K=10_000,tend=100_000,nruns=3,
        abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination"
        )

    if abs_path[1:6] == "/home/"
        abs_path = mkpath(abs_path * "/K=$K,dni=$dni,N=$N")
    else
        abs_path = abs_path * "/K=$K,dni=$dni,N=$N"
    end
    filename = "/" * filename

    #setup healthy population only initially
    n0 = healthy_pop(K)
    parameter = WrightFisher.setupparameter(model_parameter(dni,N,K), n0, tend+1)

    f = Figure(resolution=(600,300*nruns))

    for t in 1:nruns
        time = ((t-1)*tend):(t*tend)
        #setup empty population history every time
        population_history = WrightFisher.setup_pop_hist(parameter,n0,tend+1)
        #execute simulation
        WrightFisher.mainiteration!(
            n0,
            1:tend+1,
            parameter,
            population_history,
        )
        #Add cumultative statistics
        add_stats!(population_history,n0)
        #Save Data
        save(abs_path * filename *"T$t.jld",convertforsaving(population_history,N,dni,K))
        #add plot
        add_mlp_plot!(
            f,population_history,tend;
            i=t,xtimes=time,
            minload=0,minprev=0
            )
    end

    #save plot
    save(abs_path * filename *"_Toverview.pdf", f)
end

#---
