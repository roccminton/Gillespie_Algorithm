
#---
using JLD
using CairoMakie
using LinearAlgebra

#---

#calculates the forward differences between two neighbouring values of v
#and safes the result in r
function ∇!(r,v)
    n=length(v)
    !isequal(length(r),n-1) && error("Error:Dimension mismatch!")
    for (i,x) in enumerate(view(v,1:n-1))
        r[i] = v[i+1]-x
    end
end
function ∇(v)
    r = Vector{Float64}(undef,length(v)-1)
    ∇!(r,v)
    return r
end

#calculates the second order forward differences of v
#and safes the result in r
function Δ!(r,v)
    n=length(v)
    !isequal(length(r),n-2) && error("Error:Dimension mismatch!")
    for (i,x) in enumerate(view(v,1:n-2))
        r[i] = v[i+2]-2*v[i+1]+x
    end
end
function Δ(v)
    r = Vector{Float64}(undef,length(v)-2)
    Δ!(r,v)
    return r
end

#smoothes the input vector `v` by a moving mean of length `2*d` and saves the result
#in `r`
function mm!(r,v,d)
    n = length(v)
    !isequal(n-2*d,length(r)) && error("Error:Dimension Mismatch")
    for i in d+1:n-d
        r[i-d] = mean(view(v,i-d:i+d))
    end
end
function mm(v,d)
    r = Vector{Float64}(undef,length(v)-2*d)
    mm!(r,v,d)
    return r
end

function findinc(v;δ=5*10^2,ε=10.0)
    #calculate the first order forward difference, smooth is with the moving
    #mean and look for the location of the maximum value
    r = argmax(mm(∇(v),δ))
    m₁ = mean(view(v,1:r))
    m₂ = mean(view(v,r+1:length(v)))
    #if the rise is at the beginning and the means befor and after the rise
    #are close, then it was only the initla rise
    r < 100 && m₂-m₁ < ε && return 0
    return r
end

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)
prev(h) = replace_NaN(h["ill"] ./ h["PopSize"])
ml(h) = replace_NaN(h["ML"] ./ h["PopSize"])

#mean over time
mot(data, time) = vec(mean(data[:, time], dims = 2))

function findextime(h::Dict,key="PopSize")
    ext = findfirst(iszero,h[key])
    isnothing(ext) && return h["historylength"]
    return ext
end

function timestoint(ext,tend)
    if !isempty(ext)
        times = [1:ext[1]]
        for i in 2:length(ext)-1
            push!(times,(ext[i]+1):ext[i+1])
        end
        push!(times,ext[end]:tend)
    else
        times = [1:tend]
    end
    return times
end

extime(data, loadclass, birthtime, tend) = findfirst(
    iszero,
    view(data["LoadHist"], loadclass + 1, birthtime:tend),
)

extime(data::Vector, loadclass, birthtime, tend) = findfirst(
    x->iszero(x[loadclass+1]),
    view(data,birthtime:tend)
)

get_historylength(data::Dict) = data["historylength"]
get_historylength(data::Vector) = length(data)

function findextimes(data, tend = 0)
    iszero(tend) && (tend = get_historylength(data))
    extimes = [extime(data, 0, 2, tend)]
    lc = 1
    while !isnothing(extimes[end])
        ext = extime(data, lc, extimes[end], tend)
        isnothing(ext) && return extimes
        push!(extimes, ext + extimes[end])
        lc += 1
    end
end

function findcorrelationcluster(corrmx,ε=0.1)
    N = size(corrmx)[1]
    C = []
    for i in 1:N-1
        for j in i+1:N-1
            if corrmx[i,j] ≥ 1 - ε
                ci = findindinset(C,i)
                iszero(ci) ? push!(C,[i,j]) : (!(j ∈ C[ci]) && push!(C[ci],j))
            end
        end
    end
    return C
end

cluster_index(s) = iszero(length(s)) ? 0 : sum(s)
n_cluster(s) = iszero(length(s)) ? 0 : size(s)[2]
cluster_count(s) = iszero(length(s)) ? 0 : sum(count(!iszero,s,dims=1))

function CovToCor(M)
    iszero(M) && return M
    D = Diagonal(sqrt.(diag(M)))
    try
        DInv = inv(D)
        return DInv * M * DInv
    catch
        return M
    end
end

#---

function add_mlp_plot!(
    f,h,tend;
    i=1,j=1,tstart=0,
    extimes=[],dagger=[],
    maxload=0,maxprev=0,
    minload=nothing,minprev=nothing,
    align="center",xtimes=[],
    c1=:orange,c2=:red
    )

    isempty(xtimes) && (xtimes=tstart:tend)

    #create two axis which share a x-axis
    axprev = Axis(
        f[i, j],
        yticklabelcolor = c1,
        yaxisposition = :right,
        xaxisposition = :top,
        xticksvisible = false,
        ytickformat = vs -> map(v->L"$%$(round(Int,v*100))$ %",vs),
        ylabel= L"$ $Prevalence"
    )

    axload = Axis(
        f[i, j],
        yticklabelcolor = c2,
        ylabel = L"$ $Mutation Burden",
        xlabel = L"$ $Time",
        ytickformat = vs -> map(v->L"$%$(round(Int,v))$ ",vs),
        xtickformat = vs -> map(v->L"$%$(round(Int,v))$",vs)
    )

    if align == "right"
        axprev.ylabel = ""
        hideydecorations!(axprev,grid=false)
    elseif align == "left"
        axload.ylabel = ""
        hideydecorations!(axload,grid=false)
    end

    # axload.xticks = (
    #     range(xtimes[1], xtimes[end]; length = 5),
    #     [
    #         string.(round.(Integer, x / 1000)) .* "K" for
    #         x in range(xtimes[1], xtimes[end]; length = 5)
    #     ],
    # )

    if isempty(extimes)
        hidespines!(axprev)
        hidexdecorations!(axprev)
    else
        if isempty(dagger)
            axprev.xticks = (extimes,fill(L"\dagger",length(extimes)))
        else
            axprev.xticks = (extimes,dagger)
        end
    end


    #set the xlims
    linkxaxes!(axload, axprev)
    xlims!(axload, (xtimes[1], xtimes[end]))
    xlims!(axprev, (xtimes[1], xtimes[end]))

    #set the ylims
    !iszero(maxload) && ylims!(axload,high=maxload)
    !iszero(maxprev) && ylims!(axprev,high=maxprev)

    !isnothing(minload) && ylims!(axload,low=minload)
    !isnothing(minprev) && ylims!(axprev,low=minprev)

    #plot the data
    prevline = lines!(
        axprev,
        xtimes,
        prev(h)[tstart+1:tend+1],
        color = c1,
        label = "Prevalence",
    )
    loadline =
        lines!(axload, xtimes, ml(h)[tstart+1:tend+1], color = c2, label = "Mutation Load")

    return axload, axprev, loadline, prevline
end

function add_corrcluster_plot!(f,d,tend;i=1,j=1,ts=0,Δ=10)
    T = round.(Integer,tstart:Δ:tend)

    axcount = Axis(f[i,j],yaxisposition =:right,xticksvisible = false,yticklabelcolor = :green,ylabel="Genes in Cluster")
    axn = Axis(f[i,j],yticklabelcolor = :blue,ylabel="Number of Clusters",xlabel="Time")

    hidexdecorations!(axcount)

    countline = lines!(axcount,T,cluster_count.(d["CorrCluster"][T.+1]),color=:green)
    nline = lines!(axn,T,n_cluster.(d["CorrCluster"][T.+1]),color=:blue)

    xlims!(axcount,(tstart,tend))
    xlims!(axn,(tstart,tend))
    ylims!(axcount,low=0)
    ylims!(axn,low=0)

    axn.xticks = (
        range(tstart, tend; length = 5),
        [string.(round.(Integer, x / 1000)) .* "K" for x in range(0,tend;length=5)]
         )

    return axcount, anx, countline, nline
end
