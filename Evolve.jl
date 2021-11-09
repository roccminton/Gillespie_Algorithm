include("BirthDeath.jl")
include("MainFunctions.jl")

import .Gillespie
import .BirthDeath


#Setup Model parameters
b = 1.0
d = 0.9
K = 100.0
c = 10^(-5)
t = 0:150
n_0 = 10.0
model_parameter = (
        birth = b,
        death = d,
        competition = c,
        K = K
        )
config = BirthDeath.ModelConfiguration("Logistic",rescaled=true)

#execute the simulation
history = Gillespie.run_gillespie(t,n_0,model_parameter,config)

using Plots

plot(t,history)
