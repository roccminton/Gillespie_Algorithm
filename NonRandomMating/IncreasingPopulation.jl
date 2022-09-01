include("../DiploidModel.jl")
include("../Plotting.jl")

import .DiploidModel
import .PlotFromDicts

parameter_from_K(K,b,d,dni,N) = (
        birth = b,
        death = d,
        competition = (b-d)/K,
        μ = dni,
        Nloci = N
)

for K in [
        500,1_000,
        10_000,50_000,
        100_000,500_000,
        1_000_000,1_500_000,2_000_000
        ]

        Ks = [K,K*20]
        tup = 150

        b = 1.0
        d = 0.9
        dni = 1.2 * 10^(-2)
        N = 100

        parameters = [parameter_from_K(K,b,d,dni,N) for K in Ks]

        t = 0:1000

        x0 = DiploidModel.generatehealthypopulation(Ks[1],N)
        #execute the simulation
        history = DiploidModel2.runseveralgillespie(t,tup,x0,parameters)

        #plot simulation
        PlotFromDicts.plotmutationloadandprevalence(history)

        #savefig
        savefig("/Users/roccminton/Documents/Uni/Gillespie_Algorithm/PlotOutput/random_mating_$(K[1]).pdf")
end
