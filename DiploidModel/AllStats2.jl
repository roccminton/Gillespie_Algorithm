#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    loadhist :: Dict
    loadpos :: Dict
    loadpos2 :: Array
    haploidloadpos :: Vector
    covmatrix :: Vector
    corrcluster :: Vector
    clustercount :: Dict
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadpos = Dict(
            "Healthy"=>[emptytraits(par.Nloci,Int64) for _ in 1:l],
            "Ill"=>[emptytraits(par.Nloci,Int64) for _ in 1:l]
            )
    loadhist = Dict(
            "Healthy"=>[emptyhistorgram(par) for _ in 1:l],
            "Ill"=>[emptyhistorgram(par) for _ in 1:l]
            )
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    haploidloadpos = [zeros(par.Nloci+1) for _ in 1:l]
    loadpos2 = zeros(Integer,(3,par.Nloci,l))
    covmatrix = Vector{Matrix{Float64}}(undef,l)
    corrcluster = Vector{Vector{Vector{Integer}}}(undef,l)
    clustercount = Dict(
        "Healthy" => Vector{Vector{Integer}}(undef,l),
        "Ill" => Vector{Vector{Integer}}(undef,l)
        )
    return MLPLoadHistLoadPos(
        mlp,loadhist,loadpos,loadpos2,haploidloadpos,covmatrix,corrcluster,clustercount,par
        )
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
#in the main function one only knows the scope of DiploidModel2 and not Gillespie, therefore
#one needs to change and hand in the function in this module
DiploidModel2.choosestatsfunction(population_history::MLPLoadHistLoadPos) = saveonestep!

#and define your own saveonestep! function
function saveonestep!(ph::MLPLoadHistLoadPos,index,ps,par)
    #save regular MLP History
    DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)
    #save the histogram History
    savehistdata!(ph.loadpos,index,ps,par.cloadpos)
    #save the histogram History
    savehistdata!(ph.loadhist,index,ps,par.cloadhist)
    #safe the haploid histogramm data
    savelistdata!(ph.haploidloadpos,index,ps,par.chaploidloadhist)
    #save the histogram History
    savehistdata2!(ph.loadpos2,index,ps,par)
    #calculate and save covariance matrix
    savecovmatrix!(ph.covmatrix,index,ps,par,ph.loadpos)
    #calculate and save correlation cluster
    savecorrcluster!(ph.corrcluster,ph.clustercount,ph.covmatrix,index,ps,par)
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],-1,par.Nloci)
    update_haploidloadhist!(par.chaploidloadhist,par.traits[index],-1,par)
    update_histogram!(par.cloadposhist,par.traits[index],-1,par)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1,par.Nloci)
    update_loadhist!(par.cloadhist,par.traits[index],+1,par.Nloci)
    update_haploidloadhist!(par.chaploidloadhist,par.traits[index],+1,par)
    update_histogram!(par.cloadposhist,par.traits[index],+1,par)
end

#change data for better storage
function convertforsaving(h)
    #list all the parameters worth saving
    safe_parameter=[
        "death","μ","Nloci","historylength","ccuts","recombination",
        "competition", "birth", "rates", "K",
        ]
    safe_h = Dict{String,Any}()
    #safe mlp as it is
    merge!(safe_h,h.mlp)
    #loadhist
    for (key,value) in h.loadhist
        levelhists!(value)
        safe_h["LoadHist" * key] = hcat(value...)
    end
    #loadpos
    for (key,value) in h.loadpos
        safe_h["LoadPos" * key] = hcat([vcat(d...) for d in value]...)
    end
    #clustercount
    for (key,value) in h.clustercount
        levelhists!(value)
        safe_h["ClusterCount" * key] = hcat(value...)
    end
    #haploid loadhist
    safe_h["HaploidLoadHist"] = hcat(h.haploidloadpos...)
    safe_h["CovMatrixList"] = h.covmatrix
    safe_h["LoadPos2"] = h.loadpos2
    for clist in h.corrcluster
        levelhists!(clist)
    end
    safe_h["CorrCluster"] = [hcat(clist...) for clist in h.corrcluster]
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

function levelhists!(histlist)
    if !isempty(histlist)
        maxlength = maximum(length.(histlist))
        for hist in histlist
            while length(hist)<maxlength
                push!(hist,zero(eltype(hist)))
            end
        end
    end
end


#---

function DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l)
    samplesize = min(1000,par.K)
    return (
        par...,
        cloadpos = initialloadpos(par,n0),
        cloadhist = initialloadhist(par,n0),
        cloadposhist = initialhistogram(par,n0),
        chaploidloadhist = initialhaploidloadhist(par,n0),
        samplesize = samplesize,
        samplegamets = Vector{Int64}(undef,samplesize),
        sampleinds = Vector{Int64}(undef,samplesize),
        loadmeans = zeros(Float64,par.Nloci),
        DInv = zeros(Float64,(par.Nloci,par.Nloci)),
        ε = (sqrt(2)-1)/(2*sqrt(2)),
        f = fill(true,2)
        )
end

function initialloadpos(par,n0)
    hist = Dict("Healthy" => emptytraits(par.Nloci,Int64), "Ill" => emptytraits(par.Nloci,Int64))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(hist,ind,1,par.Nloci)
    end
    return hist
end

function initialloadhist(par,n0)
    hist = Dict("Healthy" => emptyhistorgram(par), "Ill" => emptyhistorgram(par))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadhist!(hist,ind,1,par.Nloci)
    end
    return hist
end

function initialhaploidloadhist(par,n0)
    hist = zeros(par.Nloci+1)
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_haploidloadhist!(hist,ind,1,par)
    end
    return hist
end

function initialhistogram(par,n0)
    hist = zeros(Integer,(3,par.Nloci))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_histogram!(hist,ind,1,par)
    end
    return hist
end

function savehistdata!(hhist,index,n0,chist)
    savelistdata!(hhist["Healthy"],index,n0,chist["Healthy"])
    savelistdata!(hhist["Ill"],index,n0,chist["Ill"])
end

function savehistdata2!(hhist,index,n0,par)
        view(hhist,:,:,index) .= par.cloadposhist
end

function savelistdata!(hlist,index,n0,clist)
    while length(hlist[index]) < length(clist)
        push!(hlist[index],zero(eltype(hlist[index])))
    end
    hlist[index] .= clist
end

function savecovmatrix!(covmatrixvec,index,n0,par,loadpos)
    #calculate means per position
    for i in 1:par.Nloci
        par.loadmeans[i] = sum(loadpos["Healthy"][index][j][i] + loadpos["Ill"][index][j][i] for j in par.choosecopyfrom)
    end
    par.loadmeans .= par.loadmeans ./ n0["PopSize"]
    #choose individuals at random
    rand!(par.sampleinds,vcat(par.indices["healthy"],par.indices["ill"]))
    #choose gamets at random
    rand!(par.samplegamets,par.choosecopyfrom)
    #calculate covariances
    C = sum(AA′(par.traits[ind_index][i_choice].-par.loadmeans) for (ind_index,i_choice) in zip(par.sampleinds,par.samplegamets))
    #normalize
    C .= C ./ (n0["PopSize"]-1)
    try
        #calculate inverse diagonal
        par.DInv .= inv(Diagonal(sqrt.(diag(C))))
        #calculate correlation matrix
        C .= par.DInv * C * par.DInv
    catch
    end
    #safe
    covmatrixvec[index] = C
end

function savecorrcluster!(corrclustervec,clustercountdict,covmatrixvec,index,n0,par)
    clusters = findcorrelationcluster(covmatrixvec[index],par.ε)
    corrclustervec[index] = clusters
    if isempty(clusters)
        clustercountdict["Healthy"][index] = fill(2*(n0["PopSize"]-n0["Ill"]),1)
        clustercountdict["Ill"][index] = fill(2*n0["Ill"],1)
    else
        #iterate over all individuals to find cluster members
        for (a,A) in [("healthy","Healthy"),("ill","Ill")]
            #setup empty count list
            clustercountdict[A][index] = zeros(Integer,length(clusters)+1)
            for ind_index in par.indices[a]
                par.f .= [true,true]
                for (n,cluster) in enumerate(clusters)
                    for i in par.choosecopyfrom
                        if par.f[i] && sum(par.traits[ind_index][i][j] for j in cluster) == length(cluster)
                             clustercountdict[A][index][n] += 1
                             par.f[i] = false
                        end
                    end
                end
                clustercountdict[A][index][end] += sum(par.f)
            end
        end
    end
end
# calculation of A*A' for a Matrix A
AA′(A) = A*A'

function update_loadpos!(hist,ind,i,Nloci)
    if DiploidModel2.ispropagable(ind,Nloci)
        hist["Healthy"] .+= i .* ind
    else
        hist["Ill"] .+= i .* ind
    end
end

function update_loadhist!(hist,ind,i,Nloci)
    load = round(Integer,DiploidModel2.mutationload(ind)+1)
    if DiploidModel2.ispropagable(ind,Nloci)
        do_update_loadhist!(hist,"Healthy",load,i)
    else
        do_update_loadhist!(hist,"Ill",load,i)
    end
end

function do_update_loadhist!(hist,key,load,i)
    while !checkbounds(Bool, hist[key],load)
        push!(hist[key],1)
    end
    hist[key][load] += i
end

function update_haploidloadhist!(hist,ind,i,par)
    for j in par.choosecopyfrom #1:2
        hist[sum(ind[j])+1] += i
    end
end

function update_histogram!(hist,ind,i,par)
    for n in 1:par.Nloci
        hist[diploid_mutation(ind,n)+1,n] += i
    end
end

#---

maxmutationload(Nloci,μ,K) = min(Nloci + quantile(Poisson(μ),1-1/K),2*Nloci)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = zeros(Int64,maxmutationload(par))

diploid_mutation(ind::Vector,n) = ind[1][n] + ind[2][n]
diploid_mutation(ind::SparseVector,n) = ind[n]

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

function findindinset(S,i)
    for (n,s) in enumerate(S)
        i ∈ s && return n
    end
    return 0
end
