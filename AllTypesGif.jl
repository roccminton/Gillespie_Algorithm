using CSV
using DataFrames
using Plots
using JLD
using SparseArrays
using Distributions
using ProgressMeter
using LaTeXStrings

include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")

import .ToolkitNRM
import .PlotFromDicts

preveq(dni) = 1-exp(-dni/2)
eqload(dni,N) = 2N*sqrt(1-exp(-dni/2N))

function dosimulation(N,dni,filename;K=10_000,tend=100_000,nruns=3)
    abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N"
    abs_path = mkpath(abs_path)
    filename = "/" * filename

    println("Currently at dni = $(round(dni,digits=2)), N = $N")
    plots = Plots.Plot[]
    for i in 1:nruns
        #Run Simulation
        h = ToolkitNRM.execute_cont(K,dni,N,tend)
        #Save Data
        save(abs_path * filename *"$i.jld",ToolkitNRM.convertforsaving(h))
        #Create Plot
        p = PlotFromDicts.plot_MLP(h.mlp)
        push!(plots,p)
    end
    #save data and plot
    savefig(plot(plots...,layout=(nruns,1),size=(600,300*nruns)),abs_path * filename * "_overview.pdf")
end

function get_points(data,rand_r,rand_θ,centers,trig_func)
    Is = vcat([0],cumsum(data))
    points = Vector{Float64}(undef,Is[end])
    for i in 1:length(Is)-1
        Is[i]≠Is[i+1] && set_points!(points,Is[i]+1:Is[i+1],centers[i],rand_r,rand_θ,trig_func)
    end
    return points
end

function set_points!(points,rd,c,rs,θs,trig_func)
    view(points,rd) .= c .+ view(rs,rd) .* trig_func.(view(θs,rd))
end

function circleshape(cx,cy,r)
    θ = LinRange(0,2π,500)
    cx .+ r*sin.(θ), cy .+ r*cos.(θ)
end

function scatter_population(data,xcenter,ycenter,xyunique=unique!(collect(zip(xcenter,ycenter))),R=1,pop_size=sum(data))
    #declare local variables
    lw = 0.15 #linewith
    fs = 20 #fontsize
    αc = 1 #alpha for circle
    αm = 0.75 #alpha for marker
    ms = 4 #markersize

    upper_textaA = ["A A", "a A", "A a", "a A", "a a", "a A", "A a", "a a", "a a", "a a"]
    lower_textaA = ["A A", "A A", "A A", "a A", "A A", "A a", "A a", "a A", "A a", "a a"]
    circlecolor = [:orange,:orange,:orange,:red,:orange,:orange,:red,:red,:red,:red]

    p = plot(aspect_ratio=1,framestyle=:zerolines,legend=false)
    for (i,(x,y)) in enumerate(xyunique)
        plot!(p,circleshape(x,y,R),color=circlecolor[i],alpha=αc)
        annotate!(p,x,y+lw,text(upper_textaA[i],:black,:center,fs))
        annotate!(p,x,y-lw,text(lower_textaA[i],:black,:center,fs))
    end

    #generate random positions for points within circle
    rs = R .* sqrt.(rand(pop_size))
    θs = 2π .* rand(pop_size)

    scatter!(p,
        get_points(data,rs,θs,xcenter,cos),
        get_points(data,rs,θs,ycenter,sin),
        aspect_ratio=1, framestyle =:zerolines,
        alpha = αm, legend = false,
        markersize = ms, markerstrokewidth = 0,
        size = (1000,1000),
        axis = false, ticks=false
        )

    return p
end

function create_gif_N2(data,abs_path,filename;giflength=1000,tend=100000)
    p = Progress(giflength)

    x_center=[0,2,-2,-1,2,3,1,2,-2,1,-3,-2,-1,2,-2,0]
    y_center=[4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4]
    xyunique=[(0,4),(-2,2),(2,2),(-3,0),(-1,0),(1,0),(3,0),(-2,-2),(2,-2),(0,-4)]

    anim = @animate for t ∈ round.(Int64,range(1,tend;length=giflength))
        #mlp
        mlp = PlotFromDicts.plot_MLP(data)
        vline!(mlp,[t],label="")
        #mutation class histogram
        hist = plot(
            xticks=(collect(1:5),collect(0:4)),ylim=(0,1),xlim=(0,6),
            xlabel="Mutation Load Classes",ylabel="Frequency")
        bar!(hist,(data["LoadHistHealthy"][:,t] .+ data["LoadHistIll"][:,t]) ./ (data["PopSize"][t]),label="healthy",color=:orange)
        bar!(hist,Vector(data["LoadHistIll"][:,t]) ./ (data["PopSize"][t]),label="ill",color=:red)
        #scatter population
        scatter_pop = scatter_population(data["Types"][:,t],x_center,y_center,xyunique)
        #layout
        l = @layout [
            [a
            b] c
        ]
        plot(mlp,hist,scatter_pop,layout=l,size=(1600,1000))
        next!(p)
    end

    return gif(anim,abs_path * "$(i)_pop.gif")
end


# K = 10000
# dni=10
# N=2
# i=1
# tend = 25000

#abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N/"
#filename = "allstats"

#h = ToolkitNRM.execute_cont(K,dni,N,tend)
#save(abs_path * filename *"_data_$i.jld",ToolkitNRM.convertforsaving(h))
#PlotFromDicts.plot_MLP(h.mlp)
#d = load(abs_path * filename *"_data_$i.jld")
#data = load(abs_path*filename*"_data_$i.jld")
#create_covmatrixgif(data["covmatrix"],data["mlp"],length(data["mlp"]["PopSize"]),abs_path,filename*"_covmatrix")


#create_gif_N2(d,abs_path,filename;giflength=10,tend=tend)
