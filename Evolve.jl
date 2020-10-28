include("BirthDeath.jl")
using .BirthDeath

include("MainFunctions.jl")

#Setup Model parameters
b = 1.0
d = 0.9
c = (b-d)/10^4
t = 0:500
n₀ = PopulationState(500)
model_parameter = Parameter(b,d,c)

#execute the simulation
#history = run_gillespie(t,n₀,model_parameter)

using Plots

plot(t,[ps.size for ps in history])
