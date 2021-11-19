include("../Multitype.jl")

import .Multitype

using LinearAlgebra
using Plots

#Setup Model parameters
ntypes = 5

#birht rates get higher for higher indexed types
b = collect(range(0.1,stop=0.12,length=ntypes))
#deaht rates stay constant
d = fill(0.09,ntypes)
#uniform competition, but only between smaller types
c = UpperTriangular(fill(10^(-4),(ntypes,ntypes)))

#Mutation from n->n+1 only
M = fill(1.0,(ntypes,ntypes))
mutation = Bidiagonal(M,:L)-Diagonal(M)
mutation[end,end] = 1.0

#time horizon
t = 0:500

#start in a monomorphic population with only the first trait
n_0 = zeros(ntypes)
n_0[1] = 100.0

#collect all model parameters
model_parameter = (
        birth = b,
        death = d,
        competition = c,
        M = mutation ,
        μ = 10^(-3)
        )

#execute the simulation
history = Multitype.rungillespie(t,n_0,model_parameter)
#execute a rescaled version of the simulation
history_res10 = Multitype.rungillespie(
        t,n_0,
        (model_parameter...,K=100.0),
        rescaled=true
        )
#execute a rescaled version of the simulation
history_res100 = Multitype.rungillespie(
        t,n_0,
        (model_parameter...,K=1000.0),
        rescaled=true
        )

#plot the simulation result
plothistory(history,ntypes) = plot(t,[history[:,i] for i in 1:ntypes],legend=false)

p = plothistory(history,ntypes)
p_res10 = plothistory(history_res10,ntypes)
p_res100 = plothistory(history_res100,ntypes)

plot(p,p_res10,p_res100,layout=(3,1))
