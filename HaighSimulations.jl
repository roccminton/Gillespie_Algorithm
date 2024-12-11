#Reproduce Haigh Simulations for haploid Muellers Ratchet

using CairoMakie
using Distributions
using ProgressMeter
using Statistics

#Gives the relative frequency of populations where untile time t J fittest populations did die out for Pₙ(τᴶ>t)
extinction_prob(Fit,J,t,n_runs) = (n_runs - Fit[t,J+1])/n_runs

function mean_extinction_time(Fit,J,n_runs)
	#if for all runs this class did die out
	if Fit[end-1,J+1] == 0
		exs = [Fit[t,J+1]-Fit[t+1,J+1] for t in 1:size(Fit)[1]-1]
		return sum( (t-1) * p for (t,p) in enumerate(exs))/n_runs
	else
		return Inf
	end
end

p(k,X,s,Poi::Distribution,T) = sum(
	X[k+2-j]*(1-s)^(k+X[1]-j) * pdf(Poi,j-X[1]) / T
	for j in X[1]:k
		)
p(k,X,s,Poi::Distribution) = p(k,X,s,Poi::Distribution,T(1,X,s))
p(k,X,s,λ) = p(k,X,s,Poisson(λ))

T(r,X,s) = sum(X[i]*(1-s)^((X[1]+i-2)*r) for i in 2:length(X))

function runHaigh(X₀,s,λ,t_end;ε=10^(-8))

	K = length(X₀)-1
	N = sum(view(X₀,2:K+1))
	Xs = Matrix{Integer}(undef,t_end,K+1)
	Xs[1,:] .= X₀
	ps = zeros(K)
	Poi = Poisson(λ)

	@showprogress for t in 2:t_end
		T₁ = T(1,view(Xs,t-1,:),s)
		ps .= [p(i,view(Xs,t-1,:),s,Poi,T₁) for i in Xs[t-1,1]:Xs[t-1,1]+(K-1)]
		while ps[end] > ε
			Xs[t-1,1] += 1
			for i in 2:K
				Xs[t-1,i] = Xs[t-1,i+1]
			end
			popfirst!(ps)
			push!(ps,p(K+Xs[t-1,1]-1,view(Xs,t-1,:),s,Poi))
		end
		Xs[t,1] = Xs[t-1,1]
		Xs[t,2:K+1] .= rand(Multinomial(N,ps))
	end

	return Xs
end

function findextime(Xs,i)
	extime = findfirst(iszero,Xs[:,i])
	if extime == 1
		birthtime = findfirst(x->x>10,Xs[:,i])
		t_end = length(Xs[:,i])
		extime = findfirst(iszero,Xs[birthtime:t_end,i])
	end
	return extime
end

allextimes(X,t,slot) = [findextime(X,i) for i in 2:findfirst(!iszero,X[t,2:slot+1])]

eq(k,N,Poi::Distribution) = N*pdf(Poi,k)
eq(k,N,θ) = eq(k,N,Poisson(θ))
mload(X) = sum(number_mutations(X,i)*X[i] for i in 2:length(X))
number_mutations(X,i) = X[1]+i-2
mean_fitness(X,s) = sum(X[i]*(1-s)^number_mutations(X,i) for i in 2:length(X))/(sum(X)-X[1])

function eq_initial(N,λ,s,slots)

	Poi = Poisson(λ/s)
	p = [pdf(Poi,k) for k in 0:slots-1]
	X₀ = vcat([0],round.(Integer,N*p))

	return X₀
end

#---

N = 10^4
s = 0.14
λ = 0.67
t_end = 25_001
slots = 50

X₀ = vcat([0,N],zeros(Integer,slots-1))
#Xs = runHaigh(X₀,s,λ,t_end)

#ML = [mload(Xs[i,:]) for i in 1:t_end] ./ N
#MF = 1 .- [mean_fitness(Xs[i,:],s) for i in 1:t_end]

extimes = [findextime(Xs,i) for i in 2:8]

f = Figure(size=(800,450),fontsize=17)

# axdagger = Axis(
# 	f[1, 1],
# 	xaxisposition = :top,
# 	yticksvisible = false,
# 	xticksvisible = false
# )
# axload = Axis(
#         f[1, 1],
#         ylabel = L"$ $Mutation Burden",
#         xlabel = L"$ $Time",
#         ytickformat = vs -> map(v->L"$%$(round(Int,v))$ ",vs),
#         xtickformat = vs -> map(v->L"$%$(round(Int,v))$",vs)
#     )
#
# axdagger.xticks = (extimes,fill(L"\dagger",length(extimes)))
# hideydecorations!(axdagger)
# ylims!(axload,low=0)
# xlims!(axload,(0,t_end))
# xlims!(axdagger,(0,t_end))
#
#
# vlines!(axload,extimes,linestyle=:dashdot,color = :black)
# lines!(axload,0:t_end-1,ML,color=RGBf(0/255,80/255,162/255))
#
# save("/home/larocca/github/Diss/Defense/haigh1_N=10000,s=0.14,l=0.67.pdf",f)
#
# f

tend = 10000

ts = vcat([1],filter(x-> x<=tend,extimes),[tend])
loadclasses = [vec(mean(Xs[ts[i]:ts[i+1],:],dims=1)) .* 10^(-3) for i in 1:length(ts)-1]
colors = [
	RGBf(6/255,42/255,119/255),
	RGBf(52/255,63/255,102/255),
	RGBf(145/255,104/255,68/255),
	RGBf(119/255,124/255,51/255),
	RGBf(237/255,145/255,34/255)
	]

axloadclass = Axis(
	f[1,1],
	ylabel = L"Frequency $\times 10^{-3}$",
	xlabel = L"$ $Classes",
	ytickformat = vs -> map(v->L"$%$(round(v,digits=1))$ ",vs),
	xtickformat = vs -> map(v->L"$%$(round(Int,v))$",vs)
	)


lins = [lines!(
	axloadclass,0:slots,lc,color=colors[i],label = iszero(i-1) ? L"$t\in [0,\dagger_{%$(i)}] " : L"$t\in [\dagger_{%$(i-1)},\dagger_{%$(i)}]"
	) for (i,lc) in enumerate(loadclasses)]


xlims!(axloadclass,(0,20))
ylims!(axloadclass,low=0)

axislegend(axloadclass)

save("/home/larocca/github/Diss/Defense/haigh2_N=10000,s=0.14,l=0.67.pdf",f)

f

#---

#eq_diploid(N,μ)=2N*sqrt(1-exp(-μ/(2N)))
#sofc(N,μ,K=10000) = K*(1-sqrt(1-exp(-μ/2N)))^(N)
#m = [abs(eq_diploid(Nloci,μ) - λ/s) for Nloci in 1:100, μ in 0.01:0.01:1]
#l = [sofc(Nloci,μ) for Nloci in 1:100, μ in 0.01:0.01:1]
#heatmap(m)
