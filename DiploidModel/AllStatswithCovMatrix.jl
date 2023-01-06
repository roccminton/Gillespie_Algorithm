#Define Type for Population History
struct MLPLoadHistLoadPos
    mlp :: Dict
    loadhist :: Dict
    loadpos :: Dict
    covmatrix :: Vector
    par :: NamedTuple
end

function DiploidModel2.setup_pop_hist(par,n₀,l)
    loadpos = Dict(
            "Healthy"=>[emptytraits(par.Nloci) for _ in 1:l],
            "Ill"=>[emptytraits(par.Nloci) for _ in 1:l]
            )
    loadhist = Dict(
            "Healthy"=>[emptyhistorgram(par) for _ in 1:l],
            "Ill"=>[emptyhistorgram(par) for _ in 1:l]
            )
    mlp = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    covmatrix = Vector{Matrix{Float64}}(undef,l)
    return MLPLoadHistLoadPos(mlp,loadhist,loadpos,covmatrix,par)
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
    #calculate and save covariance matrix
    if mod(index,100)==1
        savecovmatrix!(ph.covmatrix,index,ps,par,ph.loadpos)
    else
        ph.covmatrix[index] = zeros(par.covmatrixdim)
    end
end

function DiploidModel2.updatestats_death!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],-1)
    update_loadhist!(par.cloadhist,par.traits[index],-1)
end

function DiploidModel2.updatestats_birth!(ps,par,index)
    update_loadpos!(par.cloadpos,par.traits[index],+1)
    update_loadhist!(par.cloadhist,par.traits[index],+1)
end

#---

DiploidModel2.addstatsparameter(ph::MLPLoadHistLoadPos,par,n0,l) = (
    par...,
    cloadpos = initialloadpos(par,n0),
    cloadhist = initialloadhist(par,n0),
    covmatrixdim = (2*par.Nloci,2*par.Nloci),
    loadmeans = zeros(Float64,2*par.Nloci),
    loadstds = zeros(Float64,2*par.Nloci)
    )

function initialloadpos(par,n0)
    hist = Dict("Healthy" => emptytraits(par.Nloci), "Ill" => emptytraits(par.Nloci))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadpos!(hist,ind,1)
    end
    return hist
end

function initialloadhist(par,n0)
    hist = Dict("Healthy" => emptyhistorgram(par), "Ill" => emptyhistorgram(par))
    for ind ∈ par.traits[1:n0["PopSize"]]
        update_loadhist!(hist,ind,1)
    end
    return hist
end

function savehistdata!(hhist,index,n0,chist)
        hhist["Healthy"][index] .= chist["Healthy"]
        hhist["Ill"][index] .= chist["Ill"]
end

function savecovmatrix!(covmatrixvec,index,n0,par,loadpos)
    #setup empty cov matrix
    C = zeros(Float64,par.covmatrixdim)
    #calculate means per position
    for i in 1:par.Nloci
        par.loadmeans[i] = loadpos["Healthy"][index][1][i] + loadpos["Ill"][index][1][i]
        par.loadmeans[par.Nloci+i] = loadpos["Healthy"][index][2][i] + loadpos["Ill"][index][2][i]
    end
    par.loadmeans .= par.loadmeans ./ n0["PopSize"]
    #calculate std per position (load per position is only 0 or 1, hence first and second moment are equal)
    par.loadstds .= sqrtt.(par.loadmeans .- square.(par.loadmeans))
    #calculate covariances
    for ind_index in vcat(par.indices["healthy"],par.indices["ill"])
        for i in 1:par.Nloci
            for j in i:par.Nloci
                C[i,j] += calculate_cov(par,ind_index,1,i,i,1,j,j)/(par.loadstds[i]*par.loadstds[j])
            end
            for j in par.Nloci+1:2*par.Nloci
                C[i,j] += calculate_cov(par,ind_index,1,i,i,2,j-par.Nloci,j)/(par.loadstds[i]*par.loadstds[j])
            end
        end
        for i in par.Nloci+1:2*par.Nloci
            for j in i:2*par.Nloci
                C[i,j] += calculate_cov(par,ind_index,2,i-par.Nloci,i,2,j-par.Nloci,j)/(par.loadstds[i]*par.loadstds[j])
            end
        end
    end
    #normalize
    C .= C ./ (n0["PopSize"]-1)
    #safe
    covmatrixvec[index] = C
end

calculate_cov(par,ind_index,i_choice,i_index,i,j_choice,j_index,j) = (par.traits[ind_index][i_choice][i_index]-par.loadmeans[i])*(par.traits[ind_index][j_choice][j_index]-par.loadmeans[j])

function update_loadpos!(hist,ind,i)
    if DiploidModel2.ispropagable(ind)
        hist["Healthy"] .+= i .* ind
    else
        hist["Ill"] .+= i .* ind
    end
end

square(vec) = vec.^2
sqrtt(vec) = sqrt.(vec)

function update_loadhist!(hist,ind,i)
    load = round(Integer,DiploidModel2.mutationload(ind)+1)
    if DiploidModel2.ispropagable(ind)
        update_loadhist!(hist,"Healthy",load,i)
    else
        update_loadhist!(hist,"Ill",load,i)
    end
end

function update_loadhist!(hist,key,load,i)
    while !checkbounds(Bool, hist[key],load)
        push!(hist[key],0)
    end
    hist[key][load] += i
end

#---

maxmutationload(Nloci,μ,K) = Nloci + quantile(Poisson(μ),1-1/K^2)
maxmutationload(model_parameter) = maxmutationload(model_parameter.Nloci,model_parameter.μ,model_parameter.K)

emptyhistorgram(par) = spzeros(Integer,maxmutationload(par))
