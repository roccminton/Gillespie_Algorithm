#Collection of Haploid Wright Fisher Models

using Distributions
using SparseArrays
using Random
using Statistics
using CairoMakie
using SpecialFunctions
using ProgressMeter

function infection_prob(N,m,n)
	# There must be an overlap if n+m > N
	N < n+m && return 1.0
	# There cannot be an overlap if n=0 or m=0
	iszero(m*n) && return 0.0
	# Calculate probabilities in logspace to avoid overflow
	# Calculate logarithms of factorials
    log_N_fact = logfactorial(N)
    log_N_n_fact = logfactorial(N - n)
    log_N_m_fact = logfactorial(N - m)
    log_N_n_m_fact = logfactorial(N - n - m)

    # Calculate the numerator and denominator in log space
    numerator = log_N_n_fact + log_N_m_fact
    denominator = log_N_n_m_fact + log_N_fact

    # Calculate the probability in log space
    log_probability = numerator - denominator

    # Convert back from log space and return the result
    return 1-exp(log_probability)
end

poisson_muts(X,D::Distribution) = [sum(X[k-i+1]*pdf(D,i) for i in 0:k) for k in 0:length(X)-1]
poisson_muts(X,D::Vector) = [sum(X[k-i+1]*D[i] for i in 0:k) for k in 0:length(X)-1]
poisson_muts_scaled(X,D) = [sum(X[k-i+1]*pdf(D[k-i+1],i) for i in 0:k) for k in 0:length(X)-1]
only_one_mut_add_allways(X,p) = [X[k]*(1-p)+ (isone(k) ? 0 : X[k-1]*p) for k in 1:length(X)]
only_one_mut(X,p) = [X[k]*(1-p)+ (isone(k) ? 0 : X[k-1]*p*(1-(k-1)/N)) for k in 1:length(X)]
function two_muts(X,D)
	p,q = pdf(D,0),pdf(D,1)
	N = length(X)-1
	return [
		X[k]*p +
		(k > 1 ? X[k-1]*(1-(k-1)/N)*(q+((1-p-q)/2)*(k/N+(k-1)/N)) : 0) +
		(k > 2 ? X[k-2]*(1-p-q)*(1-(k-2)/N)*(1-(k-1)/N) : 0)
	for k in 1:length(X)
	]
end

fitness(k,X,ip) = sum(X[i]*ip[i,k] for i in 1:length(X))
fitness(X,ip) = [X[k]*fitness(k,X,ip) for k in 1:length(X)]

function generate_fitness(X,add_mut,D,ip)
	#calculate absolute fitness
	abs_fit = fitness(add_mut(X,D),ip)
	#return relative fitness
	return abs_fit ./ sum(abs_fit)
end

function run_WF(X0,add_mut,D,ip,t_end)
	K = ceil(Integer,sum(X0))
	res = [X0]
	@showprogress for t in 1:t_end
		push!(res,rand(Multinomial(K,generate_fitness(res[t],add_mut,D,ip))))
	end
	return res
end

extime(Y,k,t0) = findfirst(X->iszero(X[k]),Y[t0:end])

function extimes(Y)
	indx = 1
	t = 1
	ext = extime(Y,1,1)
	exts = Int64[]
	while !isnothing(ext)
		push!(exts,ext+t)
		t = ext
		indx += 1
		ext = extime(Y,indx,t)
	end
	return exts
end

mb(X) = sum((k-1)*X[k] for k in 1:length(X))
vb(X,mb) = sum(X[k]*(k-1 - mb)^2 for k in 1:length(X))
vb(X) = vb(X,mb(X))

#---

function runandsafe(Ns,μ;K=10_000,tend=10_000,nruns=3)
	abs_path = "/media/larocca/PortableSSD/Data/HaploidPositionFree/"
	abs_path = mkpath(abs_path)

	for N in Ns
		filename = "/N=$N,mu=$μ"
		println("Currently at N=$N")
		ip = 1 .- [infection_prob(N,n,m) for n in 0:N, m in 0:N]
		f=Figure(size=(600,300*nruns))
		for n in 1:nruns
			Y = run_WF(
				vcat([K],zeros(N)),
				poisson_muts,
				Poisson(μ),
				ip,tend
				)./K
			meanburden = mb.(Y)
			varburden = [vb(X,meanburden[i]) for (i,X) in enumerate(Y)]
	 		exs = extimes(Y)
			save(abs_path * filename * "_$n.jld", Dict(
				"HLCD" => Y, "meanMB" => meanburden,
				"varMB" =>  varburden, "extimes" => exs,
				"Nloci" => N, "μ" => μ, "K" => K, "historylength" => tend
				))
			add_mut_plot!(f,meanburden,n)
		end
	#save data and plot
	save(abs_path * filename *"_overview.pdf", f)
	end
end


function add_mut_plot!(f,mb,i)
	ax = Axis(f[i,:])
	tend = length(mb)-1

	ax.xticks = (
		range(0, tend; length = 5),
		[
			string.(round.(Integer, x / 1000)) .* "K" for
			x in range(0, tend; length = 5)
		],
	)
	xlims!(ax, (0, tend))

	lines!(ax,0:tend,mb,color=:red)
end
