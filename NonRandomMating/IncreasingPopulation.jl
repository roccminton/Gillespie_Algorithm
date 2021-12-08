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

Ks = [500,10000]
tup = 150

b = 1.0
d = 0.9
dni = 1.2 * 10^(-2)
N = 100

parameters = [parameter_from_K(K,b,d,dni,N) for K in Ks]

t = 0:1000
x0 = DiploidModel.generatehealthypopulation(Ks[1],N)

#execute the simulation
history = DiploidModel.runseveralgillespie(t,tup,x0,parameters)

#plot simulation
PlotFromDicts.plotmutationloadandprevalence(history)

#savefig
savefig("/Users/roccminton/Documents/Uni/Gillespie_Algorithm/PlotOutput/random_mating.pdf")
