include("Infinitetypes.jl")
include("MainFunctions.jl")

import .Gillespie
import .InfiniteTypes

using LinearAlgebra


#Setup Model parameters
ntypes = 10

b(x) = 1.0
d(x) = 0.9
c(x,y) = 10^(-5)

mutation = [
        0.0 0.5
        0.0 0.5
        ]
μ = 0.0
t = 0:200
n_0 = Dict(1=>[100.0,b(1),d(1)])

model_parameter = (
        birth = b,
        death = d,
        competition = c,
        compdict = InfiniteTypes.generatecompdict(n_0,c),
        M = mutation ,
        μ = μ,
        diff = 1
        )

#execute the simulation

history = Gillespie.run_gillespie(
        t,
        n_0,
        model_parameter,
        InfiniteTypes.execute!,
        InfiniteTypes.rates!,
        InfiniteTypes.setuprates(b,d,n_0),
        InfiniteTypes.setuphistory(t,n_0)
        )

using Plots

plot(t,history[1])
