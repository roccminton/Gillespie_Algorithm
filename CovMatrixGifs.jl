

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
using LinearAlgebra

include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")

import .ToolkitNRM
import .PlotFromDicts

preveq(dni) = 1-exp(-dni/2)
eqload(dni,N) = 2N*sqrt(1-exp(-dni/2N))

#generate timesteps for gif from tstart to tend of length giflength
#the total length is divided into two parts seperated at acctime with weigth accfrac for the later part
function gifsteps(tstart,tend,giflength,acctime,accfrac)
    iszero(acctime) && return round.(Integer,range(tstart,tend;length=giflength))
    l₁ = round(Integer,giflength*(1-accfrac))
    return vcat(
        round.(Integer,range(tstart,acctime;length=l₁)),
        round.(Integer,range(acctime+1,tend;length=giflength-l₁))
        )
end

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

function create_covmatrixgif(data,i,abs_path;tend=0,giflength=1000,tstart=1,filename="",acctime=0,accfrac=0.5)
    p = Progress(giflength)
    iszero(tend) && (tend=data["historylength"])
    timesteps = gifsteps(tstart,tend,giflength,acctime,accfrac)
    anim = @animate for t in timesteps
        mlp_plot = PlotFromDicts.plot_MLP(data,tend,tstart=tstart)
        vline!([t],label="",color=:green)
        covheat = heatmap(
            CovToCor(data["CovMatrixList"][t]),
            clim=(-1,1),yflip=true,xmirror=true
            )
        plot(mlp_plot,covheat,layout=(2,1),size=(600,600))
        next!(p)
    end
    return gif(anim, abs_path * filename * "covmatrix_$i.gif")
end

function CovToCor(M)
    iszero(M) && return M
    D = Diagonal(sqrt.(diag(M)))
    DInv = inv(D)
    return DInv * M * DInv
end

function create_haploidclassgif(data,i;tend=0,giflength=1000,tstart=1,filename="")

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
dni=0.6
N=60
i=1
# ε = 0.4
tend = 100_000


abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N/"
filename = "covmatrix"

#dosimulation(N,dni,filename,K=K,tend=tend,nruns=1)
#d = load(abs_path * filename *"$i.jld")
#create_covmatrixgif(d,i,abs_path)#,acctime=50000,accfrac=3/4)

# sepmutfreqs = d["LoadPosHealthy"] .+ d["LoadPosIll"]
# allmutfreqs = view(sepmutfreqs,1:d["Nloci"],:) .+ view(sepmutfreqs,d["Nloci"]+1:2*d["Nloci"],:)
# numovereps = vec(count(>=(ε*d["K"]),allmutfreqs,dims=1))
# plot(numovereps,legend=false)

t=75000
s=86000

mlp = PlotFromDicts.plot_MLP(d,tend)
vline!([44500,57500,63000,67500,75000,86000])

heatmap(
    CovToCor(mean(d["CovMatrixList"][t:s])),
    clim=(-1,1),yflip=true,xmirror=true
    )
