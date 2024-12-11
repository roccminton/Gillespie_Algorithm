using JumpProcesses
using Distributions
using CairoMakie
using DiffEqCallbacks

#=
Test implementation of example from Diekmann Doebli for TSS on continuous
trait space and with Gaussian mutation kernel
=#
struct DDParameter{F1,F2,F3,F4}
    p :: Real   #mutation probability
    m :: F1     #mutation kernel
    b :: F2     #birth rate function
    d :: F3     #death rate function
    c :: F4     #competition kernel
    K :: Real   #carrying capacity
    σ_b :: Real
    σ_c :: Real
    birthrates :: Vector #placeholder for birthrates
    deathrates :: Vector #placeholder for deathrates
    types :: Vector
end

birthrate(u, p, t) = p.birthrates[end]
function update_birthrates!(integrator,index,pm)
    diff = pm * integrator.p.b(integrator.p.types[index],integrator.p.σ_b)
    for k in index:length(integrator.p.types)
        integrator.p.birthrates[k] += diff
    end
end
function birth!(integrator)
    #choose index for birth
    i = searchsortedfirst(integrator.p.birthrates,rand()*integrator.p.birthrates[end])
    #determine mutation
    if rand() > integrator.p.p
        #in case no mutation happens at birth
        integrator.u[i] += 1
        update_birthrates!(integrator,i,1)
    else
        #in case a mutation happens at birth
        new_type = m(integrator.p.types[i])
        i = length(integrator.p.types)
        push!(integrator.p.types,new_type)
        push!(integrator.u,1)
        push!(integrator.p.birthrates,integrator.p.birthrates[end])
        push!(integrator.p.deathrates,0.0)
        integrator.u[i] += 1
        update_birthrates!(integrator,i,1)
    end
end
birth = ConstantRateJump(birthrate, birth!)

deathrate(u, p, t) = p.deathrates[end]
function death!(integrator)
    #choose index for death
    i = searchsortedfirst(integrator.p.deathrates,rand()*integrator.p.deathrates[end])
    #eliminate choosen individual
    integrator.u[i] -= 1
    #update birth rates
    update_birthrates!(integrator,i,-1)
    #update deathrates
    integrator.p.deathrates .= cumsum(
        c_t*(integrator.p.d(t) + sum(
            integrator.p.c(t,v,integrator.p.σ_c)*c_v for (v,c_v) in zip(integrator.p.types,integrator.u)
        )) for (t,c_t) in zip(integrator.p.types,integrator.u)
    )
end
death = ConstantRateJump(deathrate, death!)

#---

b(x,σ) = exp(-x^2/(2σ^2))
c(x,y,σ) = 1/1000 * exp(-(x-y)^2/(2σ^2))
d(x) = 0.0
m(x) = rand(truncated(Normal(x, 0.1), -1, 1))

par = DDParameter(
    1 / 1000^2 ,m, b, d, c,
    1000, 0.9, 0.8,
    [b(-1.0,0.9)*1000],[1000*(d(-1.0)+1000*c(-1.0,-1.0,0.8))],
    [-1.0])
u₀ = [1000]
tspan = (0.0, 1000.0)
prob = DiscreteProblem(u₀, tspan, par)
jump_prob = JumpProblem(prob,Direct(),birth,death;save_positions=(false,false))
sol = solve(jump_prob,SSAStepper(),saveat = 10.0)

#---
#plot solution
f = Figure()

ax = Axis(f[1,1])

end_index = length(sol.u)

for (i,t) in enumerate(par.types)
    start_index = findfirst(x->length(x) == i,sol.u)
    lines!(ax,
        sol.t[start_index:end],[sol.u[s][i] for s in start_index:end_index],
        color = round(Int,t*1000), colorrange = (-1000,1000)
        )
end

f
