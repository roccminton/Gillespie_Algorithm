include("../DiploidModel.jl")
include("../Plotting.jl")

import .DiploidModel
import .PlotFromDicts

function run()
        K = 1000

        b = 1.0
        d = 0.9
        c = (b-d)/(2K)

        dni = 0.1
        N = 100

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

history = run()
#plot simulation
PlotFromDicts.plotmutationloadandprevalence(history)
