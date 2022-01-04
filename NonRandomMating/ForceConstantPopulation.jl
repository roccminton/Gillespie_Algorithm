include("../DiploidModel.jl")
include("../Plotting.jl")

import .DiploidModel
import .PlotFromDicts

using Plots
using CSV

function run_constpopsize_simulation(N,dni,K,tend)

        b = 1.0
        d = 0.9

        parameter = (
                birth = b,
                death = d,
                competition = (b-d)/K,
                μ = dni,
                Nloci = N,
                rates = "allbirthrates!"
        )

        t = 0:tend
        x0 = DiploidModel.generatehealthypopulation(K,N)

        #execute the simulation
        history = DiploidModel.rungillespie(t,x0,parameter)

        abs_path = "/home/larocca/github/Gillespie_Algorithm/NonRandomMating/Output/ConstPopSize/"
        filename = "K=$K,dni=$dni,Nloci=$N"
        #convert simulation to DataFrame
        CSV.write(
                abs_path * "Data/" *filename,
                DiploidModel.historytodataframe(t,N,history)
        )

        #plot simulation
        PlotFromDicts.plotmutationloadandprevalence(history)
        #savefig
        savefig(abs_path*"Plots/"*filename*".pdf")
end

function run_forN10to100(dni)
        K = 10_000
        tend = 10_000

        for N = 10:10:100
                run_constpopsize_simulation(N,dni,K,tend)
        end
end
