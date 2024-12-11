using JLD
using StaticArrays
using ProgressMeter

include("/home/larocca/github/Gillespie_Algorithm/runandsafe_IBM.jl")
include("/home/larocca/github/Gillespie_Algorithm/ToolBoxPlotting.jl")

recs = 0.0:0.1:1.0

function run_and_mean(rec,μ,N,K,tend,b,d,c,path)
	rerun = true
	ex = nothing
	h = nothing
	H = []
	max_runs = 5
	nruns = 0
	while rerun
		h = DiploidModel2.rungillespie(
			1:tend+1,
			healthy_pop(K),
			model_parameter(K,μ,N,b,d,c,rec,"allbirthrates!",path)
		)
		#look for extinction time
		ex = findfirst(iszero,h.mlp["PopSize"])
		#if population did not goe extinct take that run
		if isnothing(ex)
			rerun = false
		#if population did go extinct but late enough take it anyway
		elseif ex > 1000
			rerun = false
		#if population did go extinct too early rerun, and save the simulation
		else
			nruns += 1
			push!(H,h)
		end
		#if the max number of runs fail without a population that survived long enough
		#break and take the average over what you got
		iszero(max_runs-nruns) && (rerun = false)
	end
	if !isnothing(ex)
		for (k,v) in h.mlp
			h.mlp[k] = v[1:ex-1]
		end
		tend = ex-1
	end
	if iszero(max_runs-nruns)
		mlp,tend = average_diff_length(H)
		return get_average(mlp,100:tend)
	else
		return get_average(h.mlp,100:tend)
	end
end

function get_average(mlp,time;ks=["PopSize","ML","Ill","C0"],fs=[1,1,1,2])
	p = mean(view(mlp[ks[1]],time))
	S = [p]
	for (i,k) in enumerate(ks[2:end])
		push!(S,mean(view(mlp[k],time))/(fs[i]*p))
	end
	push!(S,time[end])
	return S
end

function average_diff_length(H,ks=["PopSize","ML","Ill","C0"])
	ls = [length(h.mlp[ks[1]]) for h in H]
	max_length = maximum(ls)
	mlp = Dict(k => zeros(max_length) for k in ks)
	for k in ks
		for (i,h) in enumerate(H)
			view(mlp[k],1:ls[i]) .+= h.mlp[k]
		end
	end
	push!(ls,1)
	sort!(ls)
	for i in 1:length(H)
		for k in ks
			view(mlp[k],ls[i]:ls[i+1]) .*= 1/(length(H)-i+1)
		end
	end
	return mlp,max_length
end

function collect_means_varN(rec,μs,Ns,K,tend,b,d,c,abs_path,filename)
	#initialize output dictionary to save
	out = Dict(
		"rec" => rec,
		"μ" => μs,
		"N" => Ns
		)
	#initialize matrix to save means
	popsize = Matrix(undef,length(μs),length(Ns))
	ml = Matrix(undef,length(μs),length(Ns))
	prev = Matrix(undef,length(μs),length(Ns))
	#loop over all combinations
	for (i,μ) in enumerate(μs), (j,N) in enumerate(Ns)
		println("Currently at μ=$μ ($i/$(length(μs))) and N=$N ($j/$(length(Ns)))")
		res = run_and_mean(rec,μ,N,K,tend,b,d,c,abs_path*filename)
		popsize[i,j] = res[1]
		ml[i,j] = res[2]
		prev[i,j] = res[3]
	end
	#save matrices in output
	out["popsize"] = popsize
	out["ml"] = ml
	out["prev"] = prev
	#save the file
	save(abs_path * filename * "_rec=$rec.jld",out)
end

function collect_means_fixedm(rec,μ,Ns,K,tend,b,d,c,abs_path,filename)
	#initialize output dictionary to save
	out = Dict(
		"rec" => rec,
		"μ" => μ,
		"Ns" => Ns
		)
	#initialize matrix to save means
	popsize = Vector(undef,length(Ns))
	ml = Vector(undef,length(Ns))
	prev = Vector(undef,length(Ns))
	C0 = Vector(undef,length(Ns))
	extime = Vector(undef,length(Ns))
	#loop over all combinations
	for (i,N) in enumerate(Ns)
		println("Currently at μ=$μ and N=$N ($i/$(length(Ns)))")
		res = run_and_mean(rec,μ,N,K,tend,b,d,c,abs_path * filename)
		popsize[i] = res[1]
		ml[i] = res[2]
		prev[i] = res[3]
		C0[i] = res[4]
		extime[i] = res[5]
	end
	#save matrices in output
	out["popsize"] = popsize
	out["ml"] = ml
	out["prev"] = prev
	out["C0"] = C0
	out["extime"] = extime
	#save the file
	save(abs_path * filename * ".jld",out)
end

function collect_means(r,Ns,name)
	#set global variables
	K = 10_000
	tend = 10_000
	b = 1.0
	d = 0.9
	c = (b-d)/K

	abs_path = "/media/larocca/PortableSSD/Data/NoRecombination/NumericalEquilibria/"
	filename = "means_mu=0.05_rec=$(r)_" * name

	collect_means_fixedm(r,0.05,Ns,K,tend,b,d,c,abs_path,filename)
end

collect_means_smallN(r) = collect_means(r,vcat([10,15,25,50,100,150,200,250],300:100:1000),"smallNC0")
collect_means_bigN(r) = collect_means(r,[1500,2500,5000],"bigN")
collect_means_medN(r) = collect_means(r,200:100:900,"medN")
function collect_means_bigN(r,s)
	n = 10^floor(Int,log10(s))
	@showprogress dt=1 desc="Waiting..." for i in 1:n
	    sleep(s/n)
	end
	collect_means_bigN(r)
end

#d = load("/media/larocca/PortableSSD/Data/NoRecombination/SmallN/N=60,dni=0.7_2.jld")

function findloadmeans(N;μ=0.05,nruns=3)
	m_prae = 0.0
	m_post = Vector{Float64}(undef,nruns)
	extimes =  Vector{Int}(undef,nruns)
	for i in 1:nruns
		d = load("/media/larocca/PortableSSD/Data/NoRecombination/BigN/N=$N,dni=$(μ)_$i.jld")
		extimes[i] = findinc(d["ML"] ./ d["PopSize"];δ=500,ε=10)
		if iszero(extimes[i])
			m_prae += mean(d["ML"] ./ d["PopSize"])
			m_post[i] = 0
		else
			m_prae += mean(view(d["ML"],1:extimes[i]) ./ view(d["PopSize"],1:extimes[i]))
			t=max(extimes[i],d["historylength"]-10_000)
			m_post[i] = mean(view(d["ML"],t:d["historylength"]) ./ view(d["PopSize"],t:d["historylength"]))
		end
	end
	return vcat([N],extimes,m_prae / nruns,m_post)
end
