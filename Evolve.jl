include("BirthDeath.jl")
using .BirthDeath

include("MainFunctions.jl")

#Setup Model parameters
b = 1.0
d = 0.9
K = 10
c = (b-d)/(10^3)
t = [0:50,50:200]
n_0 = PopulationState(10.0)
model_parameter = [
    Parameter(
        birth = b,
        death = d,
        competition = (b-d)/500,
        K = K
        ),
    Parameter(
        birth = b,
        death = d,
        competition = (b-d)/1000,
        K = K
        )
        ]
config = fill(ModelConfiguration("Logistic",rescaled=true),2)

#execute the simulation
history = run_gillespie(t,n_0,model_parameter,config)

using Plots

plot(vcat(t...),[ps.size for ps in history])
