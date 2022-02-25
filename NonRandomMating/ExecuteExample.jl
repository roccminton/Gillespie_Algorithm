include("../DiploidModel2.jl")
include("../Plotting.jl")

import .DiploidModel2
import .PlotFromDicts
using BenchmarkTools
using Statistics
using CSV
using Plots
using DataFrames

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)

function execute_once(dni,N,K,tend)

        b = 1.0
        d = 0.9

        parameter = (
                birth = b,
                death = d,
                competition = (b-d)/K,
                μ = dni,
                Nloci = N,
                K = K,
                rates = "allbirthrates!"
        )

        t = 0:tend
        x0 = DiploidModel2.generatehealthypopulation(K)

        #execute the simulation
        history = DiploidModel2.rungillespie(t,x0,parameter)
        return history
end

function execute_and_mean(dni,N,K,tend,n)
        hs = [execute_once(dni,N,K,tend) for _ in 1:n]
        return Dict(
                k => [
                        mean(hs[i][k][t] for i ∈ 1:n)
                for t ∈ 1:tend+1]
        for k in keys(hs[1]))
end

function averageoverlast(history,key,last)
        data = history[key]
        return mean(view(data,last:length(data)))
end

Ns = 1:2
dnis = 0.01:0.01:0.2
Ks = 100:100:1000

function get_equilibrium(Ns,dnis,Ks,tend,nruns)
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
                row = [N,dni,K]
                #number of runs depend on average size
                n = ceil(Int,100_000/K)
                history=execute_and_mean(dni,N,K,tend,nruns)
                popsize = averageoverlast(history,"PopSize",last)
                push!(row,popsize)
                push!(row,averageoverlast(history,"Ill",last)/popsize)
                push!(row,averageoverlast(history,"ML",last)/popsize)
                push!(df, row)
        end
        return df
end

get_equilibrium(Ns::Integer,dnis,Ks,tend,nruns) = get_equilibrium([Ns,],dnis,Ks,tend,nruns)
get_equilibrium(Ns,dnis::Integer,Ks,tend,nruns) = get_equilibrium(Ns,[dnis,],Ks,tend,nruns)
get_equilibrium(Ns,dnis,Ks::Integer,tend,nruns) = get_equilibrium([Ns,],dnis,[Ks,],tend,nruns)

#history=execute_and_mean(dni,N,K,tend,10)
#history=execute_once(dni,N,K,tend)

#plot simulation
# PlotFromDicts.plotmutationloadandprevalence(
#         history["PopSize"],
#         replace_NaN(history["Ill"] ./ history["PopSize"]),
#         replace_NaN(history["ML"] ./ history["PopSize"]),
#         )

history = get_equilibrium(Ns,dnis,Ks,50,3)
