#=
Implementation of a simple Wright Fisher model with constant
population size
=#

using Distributions
using Random
using SpecialFunctions

stats_parameter(par,historylength) = ()

healthy_pop(K) = Dict(
	"healthy" => K,
	"ill" => 0
)

function WrightFisher.setupparameter(par, n0, historylength)

    return (
    par...,
    rndm = Vector{Int}(undef, 2),
    MutationsPerBirth = [Poisson(par.μ * (1-k/par.Nloci)) for k in 0:par.Nloci],
    infection_prob = [infection_prob(par.Nloci,n,m) for n in 0:par.Nloci, m in 0:par.Nloci],
    traits = [inittraits(par, n0),inittraits(par, n0)],
    indices = Dict(
        "healthy" => [
            collect(n0["ill"]+1:n0["ill"]+n0["healthy"]),
            collect(n0["ill"]+1:n0["ill"]+n0["healthy"]),
        ],
        "ill" => [
            collect(1:n0["ill"]),
            collect(1:n0["ill"])
        ],
    ),
    historylength = historylength,
    choosecopyfrom = 1:2,
	stats_parameter(par,historylength)...,
    )
end

function inittraits(par,n0)
    #Generate healty population with some buffer for fluctuations
    traits = [zeros(Integer,3) for _ in 1:par.K]
    #add two mutations to completely healthy individuals to get the required number of ill individuals
    for i in 1:n0["ill"]
        traits[i] .+= ones(Integer,3)
    end
    return traits
end

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

function update_stats! end

function WrightFisher.birth!(ps, par, indx, old, new)
    #choose two genetic configurations to mate
    rand!(par.rndm, par.indices["healthy"][old])
    #generate offsprings genetic configuration
    offspring!(indx, par, rand(par.MutationsPerBirth), old, new)
    #add the individual to the current population state dictionary
    if iszero(par.traits[new][indx][3])
		ps["healthy"] += 1
        #save index of new individual in appropriate list
        WrightFisher.saveorpush!(par.indices["healthy"][new], ps["healthy"], indx)
    else
        ps["ill"] += 1
        #save index of new individual in approriate list
        WrightFisher.saveorpush!(par.indices["ill"][new], ps["ill"], indx)
    end
	update_stats!(ps,par,indx,new)
end

function offspring!(offspring_index, par, n_mut, old, new)
	for i in par.choosecopyfrom
		#randomly choose one gamete and add mutation
		parental_mutations = par.traits[old][par.rndm[i]][rand(par.choosecopyfrom)]
		par.traits[new][offspring_index][i] =
			min(
				parental_mutations + rand(par.MutationsPerBirth[parental_mutations+1]),
				par.Nloci
			)
	end
	#determine if the offspring is ill or healthy
		if par.infection_prob[
			par.traits[new][offspring_index][1]+1,
			par.traits[new][offspring_index][2]+1
			] < rand()
			par.traits[new][offspring_index][3] = 0
		else
			par.traits[new][offspring_index][3] = 1
		end
	nothing
end

const_popsize(popsize,par,t) = par.K
random_fluctuations(popsize,par,t) = round(Integer,rand(par.fluctuations))


#---
