include("NonRandomMating/Toolkit.jl")
include("Plotting.jl")


import .ToolkitNRM
import .PlotFromDicts

using CSV
using DataFrames
using Plots
using Distributions
using ProgressMeter

function stringtosum(strg,start=9)
        seppos = first(findfirst(",",strg))
        return parse(Int64,strg[9:seppos-1])+ parse(Int64,strg[seppos+2:end-1])
end

function stringtovec(strg,start=9)
        seppos = first(findfirst(",",strg))
        return [parse(Int64,strg[9:seppos-1]),parse(Int64,strg[seppos+2:end-1])]
end

stringstofreq(Pos,MLP,t,N) = [stringtosum(Pos[t,n]) for n in 1:N]./MLP[t,:PopSize]
stringstofreqvec(Pos,MLP,t,N) = [stringtovec(Pos[t,n]) for n in 1:N]./MLP[t,:PopSize]

function freqtohist(freq,step=0.1,stop=0.55,start=0.05)
        hist = [count(x->x<=start,freq)]
        for ε ∈ start:step:stop-step
                push!(hist,count(x-> ε < x < ε+step, freq))
        end
        push!(hist,count(x-> x >= stop,freq))
        return hist
end

getxticks(step=0.1,stop=0.55,start=0.05) = (
        1:length(start:step:stop)+1,
        vcat(
                ["≤$(round(Integer,start*100))%"],
                ["<$(round(Integer,ε*100))%" for ε in start+step:step:stop],
                ["≥$(round(Integer,stop*100))%"]
                )
        )
#---

function creategif(dni,N,i)

        #declare local constants
        ε = 0.2
        giflength = 1000
        tend = 100000

        a = 0.05 #start for mutation per gene frequency
        b = 0.05 #step for mutation per gene frequency
        c = 0.55 #stop for mutation per gene frequency

        p = Progress(giflength) #variable for ProgressMeter

        abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=10000,dni=$dni,N=$N/"

        #load data

        MLP = DataFrame(CSV.File(abs_path * "MLP_$i.csv"))
        LoadHist_Healthy = DataFrame(CSV.File(abs_path * "LoadHist_Healthy_$i.csv"))
        LoadHist_Ill = DataFrame(CSV.File(abs_path * "LoadHist_Ill_$i.csv"))
        LoadPos_Healthy = DataFrame(CSV.File(abs_path * "LoadPos_Healthy_$i.csv"))
        LoadPos_Ill = DataFrame(CSV.File(abs_path * "LoadPos_Ill_$i.csv"))

        #data manipulation

        allfreqs = [stringstofreq(LoadPos_Healthy,MLP,t,N) .+ stringstofreq(LoadPos_Ill,MLP,t,N) for t in 1:tend]
        sepfreqs = [stringstofreqvec(LoadPos_Healthy,MLP,t,N) .+ stringstofreqvec(LoadPos_Healthy,MLP,t,N) for t in 1:tend]
        hists = freqtohist.(allfreqs,b,c,a) ./ N
        leq = filter.(x-> x < ε, allfreqs)
        geq = filter.(x-> x ≥ ε, allfreqs)
        leq_mean = mean.(leq)
        geq_mean = mean.(geq)
        leq_freq = length.(leq) ./ N
        geq_freq = length.(geq) ./ N

        xt = getxticks(b,c,a)
        hl = length(hists[1])
        popticks = floor.(Int,PlotFromDicts.eventicks(tend,4))

        #create gif
        anim = @animate for t in round.(Integer,range(1,tend;length=giflength))
                mlp = PlotFromDicts.plot_MLP(MLP)
                vline!(mlp,[t],label="")

                hist = plot(ylim=(0,0.25),xlabel="Mutation Load Classes",ylabel="Frequency")
                bar!(hist,(Vector(LoadHist_Healthy[t,:]) .+ Vector(LoadHist_Ill[t,:])) ./ (MLP[t,:PopSize]),label="healthy",color=:orange)
                bar!(hist,Vector(LoadHist_Ill[t,:]) ./ (MLP[t,:PopSize]),label="ill",color=:red)

                pos = PlotFromDicts.plot_LoadPos(sepfreqs[t],N)
                ylabel!(pos,"Loci")
                xlabel!(pos,"Frequency")

                freqhist = plot(
                        range(0,hl+1,length=tend),
                        [leq_freq,geq_freq],
                        label=["<$(round(Integer,ε*100))%" "≥$(round(Integer,ε*100))%"],
                        color= [:lightblue :salmon],
                        legend=:topleft
                        )
                vline!(freqhist,[t*(hl+1)/tend],label="")
                bar!(freqhist,
                        hists[t],
                        label="",
                        xticks=xt,rotation=25,ylabel="Percentage of genes",xlabel="Frequency of Mutation per gene",
                        ylim = (0,1),color=:blue
                )
                t′=max(1,t-1000)
                scatter!(freqhist,
                        [maximum(hists[s][i] for s ∈ t′:t) for i in 1:length(hists[1])],
                        label="", markerstrokewidth=0,color=:red
                        )

                 meanfreq = plot(
                        [leq_mean,geq_mean],
                        legend=:topleft,
                        label=["<$(round(Integer,ε*100))%" "≥$(round(Integer,ε*100))%"],
                        color= [:lightblue :salmon],
                        xlabel = "Time",
                        ylabel = "Mean Frequency per Gene",
                        xticks = (popticks,string.(popticks))
                        )
                vline!(meanfreq,[t],label="")

                l = @layout [[a; b ] c{0.3w, 0.7h} ; d e]

                plot(
                        mlp,freqhist,pos,meanfreq,hist,
                        layout=l,
                        size=(1000,800)
                        )

                next!(p)
        end

        #save and return gif
        return gif(anim, abs_path * "big_gif_$i.gif")
end
