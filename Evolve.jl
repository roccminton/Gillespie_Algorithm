include("Multitype.jl")
include("MainFunctions.jl")

import .Gillespie
import .Multitype


#Setup Model parameters
b = 1.0
d = 0.9
c = 10^(-5)
mutation = [
        0 1
        1 0
        ]
μ = 1/100
t = 0:150
n_0 = [10,100]

model_parameter = (
        birth = b,
        death = d,
        competition = c,
        M = mutation ,
        μ = μ
        )

#execute the simulation

history = Gillespie.run_gillespie(
        t,
        n_0,
        model_parameter,
        Multitype.execute!,
        Multitype.rates!
        )

using Plots

plot(t,[h[i] for h in history] for i in 1:length(n_0)])
