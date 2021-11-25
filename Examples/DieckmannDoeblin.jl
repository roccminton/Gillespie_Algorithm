include("../Infinitetypes.jl")
include("../Plotting.jl")

import .InfiniteTypes
import .PlotFromDicts

using LinearAlgebra
using Distributions

K = 1000
σ₁ = 0.9
σ₂ = 0.8

b(x, σ) = exp(-x^2 / (2σ^2))
d(x) = 0.0
c(x, y, σ, K) = inv(K) * exp(-(x - y)^2 / (2σ^2))

mutation(x) = rand(truncated(Normal(x, 0.1), -1, 1))
μ = inv(K * log(K))

model_parameter = (
        birth = x -> b(x, σ₁),
        death = d,
        competition = (x, y) -> c(x, y, σ₂, K),
        mutate = mutation,
        μ = μ,
        K = K,
)


t = 0:300
x0 = -1.0
n_x0 = InfiniteTypes.monoeq(x0, model_parameter)

#execute the simulation
history = InfiniteTypes.rungillespie(
        t,
        x0,
        n_x0,
        model_parameter,
        rescaled = true,
)
#plot simulation
PlotFromDicts.plotTSS(history)
