using GLMakie
using JumpProcesses
using Catalyst
using Distributions
using ProgressMeter

exindex(sol) = findfirst(x->iszero(x[1]),sol.u)

function mean_extime(jprob,n)
	s,real_n = 0.0, 0
	for _ in 1:n
		sol = solve(jprob,SSAStepper())
		i = exindex(sol)
		!isnothing(i) && (s += sol.t[i] ; real_n += 1)
	end
	return s/real_n
end

function median_extime!(s,jprob,n)
	for j in 1:n
		sol = solve(jprob,SSAStepper())
		i = exindex(sol)
		@inbounds(!isnothing(i) ? s[j] = sol.t[i] : s[j] = Inf)
	end
	return median(s)
end

h(N,μ) = sqrt(1-exp(-μ/N))
C(k,N,μ) = pdf(Binomial(N,h(N,μ)),k)

birthrate(u,p,t) = u[1]*p.b*exp(-p.μ/p.N)
birth!(integrator) = (integrator.u[1] += 1)
birth = ConstantRateJump(birthrate,birth!)

deathrate(u,p,t) = u[1]*(p.d + p.K * p.c)
death!(integrator) = (integrator.u[1] -= 1)
death = ConstantRateJump(deathrate,death!)

function get_av_over(μs,Ns;n_runs = 100,tspan=(0.0,1000.0),b=1.0,d=0.9,K=10_000)
	avs = Matrix{Float64}(undef,(length(μs),length(Ns)))
	u₀ = [0]
	s = Vector{Float64}(undef,n_runs)
	@showprogress for (i,μ) in enumerate(μs), (j,N) in enumerate(Ns)
		p = (b=b,d=d,c=(b-d)/K,K=K,μ=μ,N=N)
		u₀[1] = ceil(Integer,C(0,p.N,p.μ)*p.K)
		dprob = DiscreteProblem(u₀,tspan,p)
		jprob = JumpProblem(dprob,Direct(),birth,death)
		avs[i,j] = median_extime!(s,jprob,n_runs)
	end
	return avs
end

#---

μs = 0.0:0.1:1.0
Ns = 1:100

#AVs = get_av_over(μs,Ns;n_runs=100,tspan=(0.0,1000.0))

#surface(μs,Ns,AVs,axis=(type=Axis3,))

#avs = load("/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/MeanExtinctionTimes.jld")["data"]

f = Figure()

ax = Axis(f[1,1])

for i in 1:length(μs)
	lines!(ax,Ns,AVs[i,:])
end

f
