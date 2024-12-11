using CSV
using DataFrames
using Plots
using JLD
using SparseArrays
using Distributions
using ProgressMeter
using LaTeXStrings
using CSV
using Measures

include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")

import .ToolkitNRM
import .PlotFromDicts

preveq(dni) = 1-exp(-dni/2)
eqload(dni,N) = 2N*sqrt(1-exp(-dni/2N))

function dosimulation(N,dni,filename;K=10_000,tend=100_000,nruns=3,rec=0)
    abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N"
    abs_path = mkpath(abs_path)
    filename = "/" * filename

    println("Currently at dni = $(round(dni,digits=2)), N = $N")
    plots = Plots.Plot[]
    for i in 1:nruns
        #Run Simulation
        h = ToolkitNRM.execute_cont(K,dni,N,tend,rec)
        #Save Data
        save(abs_path * filename *"$i.jld",ToolkitNRM.convertforsaving(h))
        #Create Plot
        p = PlotFromDicts.plot_MLP(h.mlp)
        push!(plots,p)
    end
    #save data and plot
    savefig(plot(plots...,layout=(nruns,1),size=(600,300*nruns)),abs_path * filename * "_overview.pdf")
end

function findextimes(data,tend=0)
    iszero(tend) && (tend = data["historylength"])
    extimes = [extime(data,0,1,tend)]
    lc = 1
    while !isnothing(extimes[end])
        ext = extime(data,lc,extimes[end],tend)
        isnothing(ext) && return extimes
        push!(extimes,ext+extimes[end])
        lc += 1
    end
end

extime(data,loadclass,birthtime,tend) = findfirst(iszero,view(data["HaploidLoadHist"],loadclass+1,birthtime:tend))
birthtime(data,loadclass) = findfirst(!iszero,view(data["HaploidLoadHist"],loadclass+1,:))

function meanloadclasssize(data,loadclass,et=0)
    bt = birthtime(data,loadclass)
    isnothing(bt) && return 0
    if iszero(et)
        et = extime(data,loadclass,bt,data["historylength"])
        isnothing(et) && (et = data["historylength"])
    end
    return mean(data["HaploidLoadHist"][loadclass+1,:][bt:et])
end

C(N,μ,K,loadclass) = 2K * (1-sqrt(1-exp(-μ/N)))^(N-loadclass)
p(N,μ) = sqrt(1-exp(-μ/N))

function create_gif(data,i;tend=0,giflength=1000,tstart=1,filename="")

    p = Progress(giflength)

    iszero(tend) && (tend=data["historylength"])

    totalload = data["LoadHistHealthy"] .+ data["LoadHistIll"]
    max_freq = ceil(maximum(maximum(totalload[:,t] ./ data["PopSize"][t] for t in tstart+100:tend)),digits=2)
    haploid_max_freq = ceil(maximum(maximum(data["HaploidLoadHist"][:,t] ./ (2 .* data["PopSize"][t]) for t in tstart+100:tend)),digits=2)

    lastclass = findlast(!iszero,sum.(eachrow(totalload)))
    loadclasses = 0:lastclass-1
    lasthapclass = findlast(!iszero,sum.(eachrow(data["HaploidLoadHist"])))
    haploidloadclasses = 0:lasthapclass-1

    extimes = findextimes(data,tend)

    daggerheight = -data["K"]/10
    generange = 1:d["Nloci"]

    sepmutfreqs = data["LoadPosHealthy"] .+ data["LoadPosIll"]
    allmutfreqs = view(sepmutfreqs,generange,:) .+ view(sepmutfreqs,data["Nloci"]+1:2*data["Nloci"],:)

    mutpos_max_freq = ceil(maximum(allmutfreqs/data["K"]),digits=2)

    #create gif
    anim = @animate for t in round.(Integer,range(tstart,tend;length=giflength))
            mlp = PlotFromDicts.plot_MLP(data,tend,tstart=tstart)
            vline!(mlp,[t],label="",color=:green)
            #mutation class histogram
            hist = plot(
                ylim=(0,max_freq),
                xlabel="Mutation Load Classes",ylabel="Frequency"
                )
            bar!(hist,loadclasses,view(totalload,1:lastclass,t) ./ (data["PopSize"][t]),label="healthy",color=:orange)
            bar!(hist,loadclasses,view(data["LoadHistIll"],1:lastclass,t) ./ (data["PopSize"][t]),label="ill",color=:red)

            #haploid mutation class histogram
            haploidhist = plot(
                xlabel="Haploid Mutation Load Classes",ylabel="Frequency",legend=false,ylim=(0,haploid_max_freq))
                bar!(haploidloadclasses,view(data["HaploidLoadHist"],1:lasthapclass,t) ./ (2 .* data["PopSize"][t]),color=:blue)

            #mutation load classes
            mutpos = plot(
                xlim=(0,data["Nloci"]+1),ylim=(0,mutpos_max_freq),
                xlabel = "Loci", ylabel = "Average Mutations",
                legend = false
            )
            bar!(
                mutpos,
                generange,
                allmutfreqs[:,t]./data["PopSize"][t],
                label="",color=:blue,orientation=:v
                )
            hline!([0.4],color=:red,label="")

            #vertical bars for extinction times
            if !isnothing(extimes) && t ≥ extimes[1]
                ct = findfirst(x->x>t,extimes)
                vline!(
                    mlp,
                    isnothing(ct) ? extimes : view(extimes,1:ct-1),
                    label="",color=:black,lw=1
                )
                annotate!(
                    mlp,
                    isnothing(ct) ? extimes : view(extimes,1:ct-1),
                    daggerheight,
                    text("†",:black,:center,10)
                    )
            end

            emptyplot = plot(frame=:none)

            l = @layout [
                [a
                b{0.3h}] grid(2,1){0.3w}
            ]

            plot(
                margin = 15mm,
                mlp,mutpos,hist,haploidhist,
                layout=l,
                size=(1748,1240) #A6 quer
            )

            next!(p)
    end

    #save and return gif
    return gif(anim, abs_path * filename * "happosgif_$i.gif")
end

K = 10_000
dni=0.4
N=90
i=1
# ε = 0.4
# tend = 100_000
#
#h = ToolkitNRM.execute_cont(K,dni,N,tend)


abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N/"
filename = "corrcluster"


#dosimulation(N,dni,filename,K=K,tend=tend,nruns=1)

#d = load(abs_path * filename *"$i.jld")

#create_gif(d,i,filename="zoom_",tstart=7500,tend=27500,giflength=500)

# sepmutfreqs = d["LoadPosHealthy"] .+ d["LoadPosIll"]
# allmutfreqs = view(sepmutfreqs,1:d["Nloci"],:) .+ view(sepmutfreqs,d["Nloci"]+1:2*d["Nloci"],:)
# numovereps = vec(count(>=(ε*d["K"]),allmutfreqs,dims=1))
# plot(numovereps,legend=false)

# t = 500
#
# @gif for t in round.(Integer,range(1,tend;length=500))
#     bar([
#         (h.loadpos[1,:,t] .+ h.loadpos[2,:,t] .+ h.loadpos[3,:,t]) ./ h.mlp["PopSize"][t],
#         (h.loadpos[2,:,t] .+ h.loadpos[3,:,t]) ./ h.mlp["PopSize"][t],
#         (h.loadpos[3,:,t]) ./ h.mlp["PopSize"][t],
#         ])
# end

#----

function ClusterSizes(data)
    #maximum number of clusters
    nc = maximum(size(data["CorrCluster"][t])[end] for t in 1:data["historylength"])
    #setup empty matrix
    cs = zeros(Integer,(data["historylength"],nc))
    for (t,cc) in enumerate(data["CorrCluster"])
        !iszero(length(cc)) && (view(cs,t,1:size(cc)[2]) .= vec(reshape(count(!iszero,cc,dims=1),:,1)))
    end
    return cs
end

function reduce_data(data)
    small_d = Dict()

    #keys that just get copied
    keys = ["μ","Nloci","historylength","K","PopSize","HaploidLoadHist"]
    for key in keys
        small_d[key] = data[key]
    end

    small_d["Prev"] = PlotFromDicts.replace_NaN(d["Ill"] ./ d["PopSize"])
    small_d["ml"] = PlotFromDicts.replace_NaN(d["ML"] ./ d["PopSize"])

    small_d["ClusterCount"] = ClusterSizes(data)

    small_d["LoadPos00"] = data["LoadPos2"][1,:,:]
    small_d["LoadPos01"] = data["LoadPos2"][2,:,:]
    small_d["LoadPos11"] = data["LoadPos2"][3,:,:]

    return small_d
end

save(abs_path * filename *"$(i)_small.jld",reduce_data(d))
