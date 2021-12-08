include("../DiploidModel.jl")
include("../Plotting.jl")

import .DiploidModel
import .PlotFromDicts

function run(N,dni)
        K = 10_000

        b = 1.0
        d = 0.9
        c = (b-d)/(K)

        model_parameter = (
                birth = b,
                death = d,
                competition = c,
                μ = dni,
                Nloci = N
                )

        t = 0:100
        x0 = DiploidModel.generatehealthypopulation(K,N)

        #execute the simulation
        history = DiploidModel.rungillespie(t,x0,model_parameter)

        return history
end

N = 50

for dni ∈ 1.2:-0.1:0.1
        history = run(N,dni)
        #plot simulation
        PlotFromDicts.plotmutationloadandprevalence(history)
        #safefig
        savefig("PlotOutput/N=$N;dni=$dni.pdf")
end
