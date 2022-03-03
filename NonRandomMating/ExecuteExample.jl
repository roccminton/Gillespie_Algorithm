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
        return DiploidModel2.rungillespie(t,x0,parameter)
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

function get_equilibrium(Ns,dnis,Ks,tend)
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
                row = [N,dni,K]
                #number of runs depend on average size
                n = ceil(Int,100_000/K)
                history=execute_and_mean(dni,N,K,tend,n)
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

#----------------

Ns = vcat(1:10,15:5:50,60:10:100)
dnis = 0.1:0.1:1.2
Ks = vcat(100:100:1_000,1_500:500:10_000,11_000:1_000:50_000,55_000:5_000:100_000)

history = get_equilibrium(Ns,dnis,Ks,100)

abs_path = "/home/larocca/github/Gillespie_Algorithm/NonRandomMating/Output/ConstPopSize/"
#safe DataFrame
CSV.write(
        abs_path * "Data/" * "K=$(Ks[1])-$(Ks[end]),dni=$(dnis[1])-$(dnis[end]),N=$(Ns[1])-$(Ns[end])",
        history
)