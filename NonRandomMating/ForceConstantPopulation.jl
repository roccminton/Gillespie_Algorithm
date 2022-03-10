include("../DiploidModel2.jl")
include("../Plotting.jl")

import .DiploidModel2
import .PlotFromDicts

using Plots
using CSV

function execute_once(N,dni,K,tend)
        x0 = DiploidModel2.generatehealthypopulation(K)
        execute_once(N,dni,K,x0,tend)
end

function execute_once(N,dni,K,x0,tend)

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

        #execute the simulation
        return DiploidModel2.rungillespie(t,x0,parameter)
end

function execute_once(N,dni,K,x0,tend)

function execute_and_mean(dni,N,K,tend,n)
        hs = [execute_once(dni,N,K,tend) for _ in 1:n]
        return Dict(
                k => [
                        mean(hs[i][k][t] for i ∈ 1:n)
                for t ∈ 1:tend+1]
        for k in keys(hs[1]))
end

function save_data(abs_path,filename,history,t)
        #convert simulation to DataFrame
        CSV.write(
                abs_path * "Data/" *filename,
                DiploidModel2.historytodataframe(t,history)
        )
end

function safe_plot(abs_path,filename,history)
        #plot simulation
        plot_MLP(history)
        #savefig
        savefig(abs_path*"Plots/"*filename*".pdf")
end

filename(K,dni,N) = "K=$K,dni=$dni,Nloci=$N"

plot_MLP(history) = PlotFromDicts.plotmutationloadandprevalence(
                history["PopSize"],
                history["Ill"] ./ history["PopSize"],
                history["ML"] ./ history["PopSize"])

#abs_path = "/home/larocca/github/Gillespie_Algorithm/NonRandomMating/Output/ConstPopSize/"
abs_path = "/Users/roccminton/Documents/Uni/Gillespie_Algorithm/Output/"

#history = execute_once(10,1.0,10_000,500)
plot_MLP(history)
