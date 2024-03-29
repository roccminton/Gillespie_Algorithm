"""

    module ToolkitNRM

    Helpful functions to for executing simulations and plotting
    in the NRM-Framework for Gillespie and Wright-Fisher model.

"""

module ToolkitNRM

include("../DiploidModel2.jl")
include("../Plotting.jl")
include("../WrightFisher.jl")

import .DiploidModel2
import .PlotFromDicts
import .WrightFisher

using Statistics
using CSV
using DataFrames
using Plots
using SparseArrays

"""
Returns a named tuple parameter set for Gillespie algorithm
"""
cont_parameter(K,b,d,dni,N,rec=1) = (
        birth = b,
        death = d,
        competition = (b-d)/K,
        μ = dni,
        Nloci = N,
        K = K,
        rates = "allbirthrates!",
        recombination = rec
)

cont_parameter(K,dni,N,rec=1) = cont_parameter(K,1.0,0.9,dni,N,rec)

wf_parameter(K,dni,N) = (μ=dni, Nloci=N, K=K, popsize=WrightFisher.const_popsize)

"""
Executes one run of the Gillespies algorithm with fixed birth
and death rate and handed in initial population.
"""
execute_once(parameter,x0,time,alg) = alg(time,x0,parameter)

execute_wf(K,dni,N,tend) = execute_once(
        wf_parameter(K,dni,N),
        DiploidModel2.generatehealthypopulation(K),
        1:tend+1,
        WrightFisher.run_wrightfisher
)

execute_cont(K,dni,N,tend) = execute_once(
        cont_parameter(K,dni,N),
        DiploidModel2.generatehealthypopulation(K),
        0:tend,
        DiploidModel2.rungillespie
)

"""
Executes several runs of Gillespies Algorithm and takes the
average over them.
"""
function execute_and_mean(K,dni,N,tend,n,ex_func)
        hs = [ex_func(K,dni,N,tend) for _ in 1:n]
        return Dict(
                k => [
                        mean(hs[i][k][t] for i ∈ 1:n)
                for t ∈ 1:tend+1]
        for k in keys(hs[1]))
end

"""
Saves the handed in data generated by Gillespies Algorithm
as a csv file at the given location.
"""
function save_history(abs_path,filename,history)
        #convert simulation to DataFrame
        CSV.write(
                abs_path * "Data/" *filename,
                DiploidModel2.historytodataframe(history)
        )
end

function save_histogram(abs_path,filename,history)
        #convert simulation to DataFrame
        CSV.write(
                abs_path * "Data/" *filename,
                histogramtodataframe(history)
        )
end

function histogramtodataframe(history)
        df = DiploidModel2.historytodataframe(dictfromhist(history))
        df[!,:MutationLoadHistHealthy] = history["Healthy"]
        df[!,:MutationLoadHistIll] = history["Ill"]
        return df
end

"""
Generates the Mutation Load and Prevalence Plot from a history
file and saves it under the specific location.
"""
function safe_plot(abs_path,filename,history)
        #plot simulation
        plot_MLP(history)
        #savefig
        savefig(abs_path*"Plots/"*filename*".pdf")
end

filename(K,dni,N) = "K=$K,dni=$dni,Nloci=$N"

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)

"""
Generates the Mutation Load and Prevalence Plot.
"""
plot_MLP(history) = PlotFromDicts.plotmutationloadandprevalence(
                history["PopSize"],
                replace_NaN(history["Ill"] ./ history["PopSize"]),
                replace_NaN(history["ML"] ./ history["PopSize"]))
plot_MLP(df::DataFrame) = PlotFromDicts.plotmutationloadandprevalence(
                df.PopSize,
                replace_NaN(df.Ill ./ df.PopSize),
                replace_NaN(df.Mutation ./ df.PopSize))


"""
Generates new data averaging over the last entries for
every key in the history dictionary.
"""
function averageoverlast(history,key,last)
        data = history[key]
        return mean(view(data,last:length(data)))
end

"""
Executes a number n of simulations for every combination of Ns, dnis,
Ks up to tend and the calculates the average ofer the last tend/2 time
steps and over all simulations with the same parameter.
"""
function get_equilibrium(Ns,dnis,Ks,tend,ex_func=execute_cont)
        #create data frame
        df = DataFrame(
                Nloci = Int[],
                dni = Float64[],
                K = Int[],
                PopSize=Float64[],
                Prev=Float64[],
                MutationLoad=Float64[]
        )
        #set timespan to average over for equilibrium
        last = ceil(Int,tend/2)
        #loop over all instances
        for N ∈ Ns, dni ∈ dnis, K ∈ Ks
                println("Currently at N=$N,dni=$dni,K=$K")
                row = Float64[N,dni,K]
                #number of runs depend on average size
                n = ceil(Int,10_000/K)
                n=1
                history=execute_and_mean(K,dni,N,tend,n,ex_func)
                popsize = averageoverlast(history,"PopSize",last)
                push!(row,popsize)
                push!(row,averageoverlast(history,"Ill",last)/popsize)
                push!(row,averageoverlast(history,"ML",last)/popsize)
                push!(df, row)
        end
        return df
end

get_equilibrium(Ns::Number,dnis,Ks,tend) = get_equilibrium([Ns,],dnis,Ks,tend)
get_equilibrium(Ns::Number,dnis::Number,Ks,tend) = get_equilibrium([Ns,],[dnis,],Ks,tend)
get_equilibrium(Ns::Number,dnis,Ks::Number,tend) = get_equilibrium([Ns,],dnis,[Ks,],tend)
get_equilibrium(Ns,dnis::Number,Ks,tend) = get_equilibrium(Ns,[dnis,],Ks,tend)
get_equilibrium(Ns,dnis::Number,Ks::Number,tend) = get_equilibrium(Ns,[dnis,],[Ks,],tend)
get_equilibrium(Ns,dnis,Ks::Number,tend) = get_equilibrium(Ns,dnis,[Ks,],tend)

loadfromhist(hist, key) = [
    iszero(count(!iszero, h)) ? 0 :
    sum((load - 1) * freq for (load, freq) in zip(findnz(h)...)) for
    h in hist[key]
]

sizefromhist(hist, key) = [sum(findnz(h)[2]) for h in hist[key]]

function dictfromhist(hist)
    ill = sizefromhist(hist, "Ill")
    popsize = sizefromhist(hist,"Healthy") .+ ill
    mutation = loadfromhist(hist, "Ill") .+ loadfromhist(hist, "Healthy")
    return Dict(
    "PopSize" => popsize,
    "Ill" => ill,
    "ML" => mutation,
    "cutsat" => hist["cutsat"]
    )
end

end #module ToolkitNRM
