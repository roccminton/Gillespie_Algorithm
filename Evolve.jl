include("BirthDeath.jl")
using .BirthDeath

include("MainFunctions.jl")

#Setup Model parameters
b = 1.0
d = 0.9
K = 100.0
c = 10^(-5)
t = 0:150
n_0 = PopulationState(10.0)
model_parameter = Parameter(
        birth = b,
        death = d,
        competition = c,
        K = K
        )
config = ModelConfiguration("Logistic",rescaled=true)

#execute the simulation
history = run_gillespie(t,n_0,model_parameter,config)

using Plots

plot(t,[ps.size for ps in history])
