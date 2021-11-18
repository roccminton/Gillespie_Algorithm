include("Multitype.jl")
include("MainFunctions.jl")

import .Gillespie
import .Multitype

using LinearAlgebra


#Setup Model parameters
ntypes = 10

b = fill(1.0,ntypes)
d = fill(0.9,ntypes)
c = Diagonal(fill(10^(-5),ntypes))
mutation = [
        0.0 0.5
        0.0 0.5
        ]
μ = 0.0
t = 0:200
n_0 = fill(10,ntypes)

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

plot(t,[history[:,i] for i in 1:length(n_0)])
