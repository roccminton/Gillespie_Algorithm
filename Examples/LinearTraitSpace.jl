include("../Multitype.jl")
include("../MainFunctions.jl")

import .Gillespie
import .Multitype

using LinearAlgebra
using Plots

#Setup Model parameters
ntypes = 5

#birht rates get higher for higher indexed types
b = collect(range(1.0,stop=1.2,length=ntypes))
#deaht rates stay constant
d = fill(0.9,ntypes)
#uniform competition, but only between smaller types
c = UpperTriangular(fill(10^(-5),(ntypes,ntypes)))

#Mutation from n->n+1 only
M = fill(1.0,(ntypes,ntypes))
mutation = Bidiagonal(M,:L)-Diagonal(M)
mutation[end,end] = 1.0
#at a rate of aprrox 1/K
μ = 10^(-5)

#time horizon
t = 0:500

#start in a monomorphic population with only the first trait
n_0 = zeros(Int,ntypes)
n_0[1] = 100

#collect all model parameters
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

#plot the simulation result
plot(t,[history[:,i] for i in 1:ntypes],legend=false)
