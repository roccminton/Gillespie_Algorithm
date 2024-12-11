using JLD
using SparseArrays
using ProgressMeter
using Statistics
using CairoMakie

include("/home/larocca/github/Gillespie_Algorithm/ToolBoxPlotting.jl")

#d = load("/media/larocca/PortableSSD/Data/NoRecombination/CorrMatrix/N=600,dni=0.05/data_1.jld")
#d = load("/media/larocca/PortableSSD/Data/NoRecombination/BigN/N=600,dni=0.05_2.jld")
#C = load("/media/larocca/PortableSSD/Data/NoRecombination/CorrMatrix/N=500,dni=0.06/data_1CovMatrixs.jld")
#

prev(h) = replace_NaN(h["Ill"] ./ h["PopSize"])

function HLCD_before_after(d)
    f = Figure(size=(800,450),fontsize=17)

    s1 = 250
    e1 = findfirst(x->!in(1,x.nzind),d["HaploidLoadHist"])

    s2 = 80_000
    e2 = d["historylength"]

    hlcs(s,e) = sum(d["HaploidLoadHist"][s:e] ./ (2*d["PopSize"][s:e])) ./ (e-s)
    hlc1,hlc2 = hlcs(s1,e1),hlcs(s2,e2)
    m1 = maximum(hlc1.nzind)
    m2 = maximum(hlc2.nzind)
    m = max(m1,m2)
    ax = Axis(
        f[1,1],
        ytickformat = vs -> map(v->L"$%$(round(Int,v*100))$ %",vs),
        xtickformat = vs -> map(v->L"$%$(round(Int,v))$",vs),
        ylabel = L"$ $Gamete Fraction",
        xlabel = L"$ $Number of Mutations",
    )
    b1 = barplot!(ax,0:m-1,hlc1[1:m],label=L"$t \in T_1$",color=RGBf(6/255,42/255,119/255))
    b2 = barplot!(ax,0:m-1,hlc2[1:m],label=L"$t \in T_2 $",color=RGBf(237/255,145/255,34/255))
    ylims!(ax,low=0)
    axislegend(ax,)

    #save("/home/larocca/github/Diss/Defense/HLCD_N=500,dni=0.06.pdf", f)

    return f
end

sortcormatrix(corr,cluster) = corr[cluster,cluster]
sortcormatrix(corr) = sortcormatrix(corr,findcluster(corr))
sortloadpos(loadpos,cluster) = loadpos[cluster]
time_to_sampleinds(T,times) = searchsortedlast(times,T[1]):searchsortedfirst(times,T[end])
function meanr2!(L,P,d,N,l,Ts)
    L .= 0.0
    @showprogress map(t->addr2!(L,P,d,t,N,l),Ts)
    L .= L ./ length(Ts)
end

function meanr2(d,N,Ts)
    L = zeros(N,N)
    P = zeros(N)
    meanr2!(L,P,d,N,d["samplesize"],Ts)
    return L
end

function addr2!(L,P,d,t,N,l)
    get_Ps!(P,d,t,N,l)
    for n in 1:N
        if iszero(P[n])
            rnn = 0.0
        else
            rnn = D(P[n],P[n],P[n],norm=false)^2 / (P[n]*(1-P[n]))^2
            for m in n+1:N
                if iszero(P[m])
                    rnm = 0.0
                else
                    pnm = sum(view(d["SnapShots"],n,:,t) .* view(d["SnapShots"],m,:,t)) ./ l
                    rnm = D(P[n],P[m],pnm,norm=false)^2 / (P[n]*(1-P[n]) * P[m]*(1-P[m]))
                end
                L[n,m] += rnm
                L[m,n] += rnm
            end
        end
        L[n,n] += rnn
    end
end

function get_Ps!(P,d,t,N,l)
    for n in 1:N
        P[n] = sum(view(d["SnapShots"],n,:,t)) ./ l
    end
end

function D(pn,pm,pnm;norm=true)
    D = pnm-pn*pm
    iszero(D) && return D
    return (norm ? D/D_max(D,pn,pm) : D)
end

function findcluster(corr,ε=0.3)
    N = size(corr)[1]
    cluster = []
    for c in 1:N-1
        for r in c+1:N
            if corr[c,r] ≥ 1-ε
                i = findfirst(x->c∈x,cluster)
                if isnothing(i)
                    push!(cluster,[c,r])
                elseif r ∉ cluster[i]
                    push!(cluster[i],r)
                end
            end
        end
    end
    #push!(cluster,filter(x->x∉reduce(vcat,cluster),1:N))
    return cluster
    return reduce(vcat,cluster)
end

function AlleleFrequency(d)

    extimes = findextimes(d)
    ts = 80_000:d["historylength"]
    N = d["Nloci"]


    # maximum number of mutation per haploid genome
    # N_max = maximum(findlast(!iszero,d["LoadHist"][:,t]) for t in ts)
    # #collect corrmatrix and loadpositions
    # println("Calculating covariance matrices ... ")
    # cormatrix = meanr2(d,N,time_to_sampleinds(ts,d["savesnap"]))
    # loadpos = mean(d["LoadPos"][:, ts], dims = 2) ./ d["K"]
    # #reorder matrix and vector to better see cluster
    # cluster = findcluster(cormatrix)
    # cormatrix = sortcormatrix(cormatrix,cluster)
    # loadpos = sortloadpos(loadpos,cluster)

    # cluster = sort(filter(x-> length(x) > 1, findcluster(cormatrix)),by=length)
    # ls = [length(c) for c in cluster]
    # cl = vcat([0],cumsum(length(c) for c in cluster))
    #
    f=Figure(size=(800,400))
    ax = Axis(
        f[1,1],
        ytickformat = vs -> map(v->L"$%$(round(v*100,digits=1))$ %",vs),
        xtickformat = vs -> map(v->L"$%$(round(Int,v))$",vs),
        ylabel = L"$ $ Allele Frequency",
        xlabel = L"$ $ Gene Position (reordered)"
        )

    vlines!(ax,cl,linestyle=:dash,linewidth=1,color=:black)

    cmin = minimum(filter(x-> x>1,ls))
    cmax = maximum(ls)

    for i in 1:length(cluster)
        barplot!(
            ax,
            cl[i]+1:cl[i+1],
            loadpos[cluster[i]],
            color=ls[i],
            colormap=:viridis,
            colorrange=(cmin,cmax)
            )
    end

    ylims!(ax,(0.15,0.175))
    xlims!(ax,(-1,601))

    Colorbar(
        f[:,2],limits=(cmin,cmax),colormap=:viridis,label=L"$ $ Cluster Size",
        tickformat = vs -> map(v->L"$%$(round(Int,v))$",vs)
    )

    save("/home/larocca/github/Diss/Figures/AlleleFrequency_sorted_N=600,dni=0.05.pdf", f)

    f
end

# nclus = []
# ngenes = []
# T = 1:1000:d["historylength"]-1000
# for t in T
#     c = meanr2(d,N,time_to_sampleinds(t:t+1000,d["savesnap"]))
#     cs = findcluster(c)
#     l = length(cs)
#     push!(nclus,l)
#     push!(ngenes,iszero(l) ? 0 : sum(length(c) for c in cs))
# end

f = Figure(size=(800,450))

# axload, axprev, loadline, prevline = add_mlp_plot!(
#     f,d,80_000,
#     extimes = [],
#     dagger = [],
#     c2 = RGBf(0/255,80/255,162/255),
#     c1 = RGBf(198/255,38/255,6/255)
#     )

#create two axis which share a x-axis
# axcount = Axis(
#     f[2, 1],
#     yticklabelcolor = :green,
#     yaxisposition = :right,
#     xaxisposition = :top,
#     xticksvisible = false,
#     ytickformat = vs -> map(v->L"$%$(round(Int,v))$ ",vs),
#     ylabel= L"$ $Number of Cluster"
# )
#
# axgenes = Axis(
#     f[2, 1],
#     yticklabelcolor = :blue,
#     ylabel = L"$ $Genes in Cluster",
#     xlabel = L"$ $Time",
#     ytickformat = vs -> map(v->L"$%$(round(Int,v))$ ",vs),
#     xtickformat = vs -> map(v->L"$%$(round(Int,v))$",vs)
# )
#
# hidespines!(axcount)
# hidexdecorations!(axcount)
#
# T = 1000:1000:d["historylength"]
#
# t1 = 25_000
# t2 = 31_000
# t3 = 49_500
#
# vlines!(axgenes,[t1,t2,t3],color=:black,linewidth=1,linestyle=:dash)
# vlines!(axload,[t1,t2,t3],color=:black,linewidth=1,linestyle=:dash)
#
# meancount = vcat(
#     fill(14,t2-t1+1),
#     fill(13,t3-t2),
#     fill(12,d["historylength"]-t3)
#     )
#
# lines!(axcount,t1:d["historylength"],meancount,color=:darkgreen,linewidth=1)
#
# lines!(axcount,T,round.(Int,nclus),color=:green,label="Number of Cluster")
# lines!(axgenes,T,round.(Int,ngenes),color=:blue,label="Genes in Cluster")
#
# xlims!(axcount,(0,75_000))
# xlims!(axgenes,(0,75_000))
# ylims!(axcount,low=0)
# ylims!(axgenes,low=0)
#
# p=prev(d)
#
# meanp = vcat(
#     fill(mean(p[t1:t2]),t2-t1+1),
#     fill(mean(p[t2:t3]),t3-t2),
#     fill(mean(p[t3:end]),d["historylength"]-t3)
#     )
# lines!(axprev,t1:d["historylength"],meanp,color=:darkorange,linewidth=3)
#
# save("/home/larocca/github/Diss/Figures/clustersize_N=600,dni=0.05.pdf", f)

HLCD_before_after(d)
