
generatehealthypopulation(popsize) = initpopulation(popsize,0,0)

function initpopulation(popsize,ML,Ill)
    if ML >= 2*Ill
        return Dict(
            "PopSize" => popsize,
            "Ill" => Ill,
            "ML" => ML
            )
    else
        error("Mutation Load ($ML) must be at least twice the number of ill individual($Ill)")
    end
end

function timesfromswitches(totaltime,switches)
    times = Vector{typeof(totaltime)}(undef,length(switches)+1)
    tstart = 1
    for (i,t) in enumerate(switches)
        tend = findfirst(x -> x≥t,totaltime)
        times[i] = totaltime[tstart:tend]
        tstart = tend+1
    end
    times[end] = totaltime[tstart:end]
    return times
end

#--- common functions for statistics and saving

#creates a text file with all parameters worth saving
function saveparameter(par,abs_path)
    open(abs_path,"w+") do io
        println(io, "K = $(par.K)")
        println(io, "μ = $(par.μ)")
        println(io, "Nloci = $(par.Nloci)")
        println(io, "recombination rate = $(par.recombination)")
        if par.recombination ≠ 1
            if length(par.ccuts) > 1
                cutsat = join([c[end] for c in par.ccuts[1:end-1]])
                println(io, "cutsat = $cutsat")
            else
                println(io, "cutsat = no chromosome cuts")
            end
        end
    end
end

#create path in no recombination system
function create_path(abs_path, par)
    return mkpath(abs_path * "/K=$(par.K),dni=$(par.μ),N=$(par.Nloci)")
end

#--- common functions for statistics and saving

function mlptodataframe(history,folder_abs_path,name_ext = "")
    df = DataFrame(
        PopSize=history.mlp["PopSize"],
        Ill=history.mlp["Ill"],
        Mutation=history.mlp["ML"],
        )

    CSV.write(folder_abs_path * "/MLP" * name_ext * ".csv", df)

    return df
end

function loadhisttodataframe(history, folder_abs_path, name_ext = "")
    max_load = PlotFromDicts.maxmutatioins(history.loadhist)
    tend = history.par.historylength
    df_ill = empty_dataframe(max_load,tend,zero(Integer),0)
    df_healthy = empty_dataframe(max_load,tend,zero(Integer),0)

    populate_df_loadhist!(df_ill,history.loadhist["Ill"],tend)
    populate_df_loadhist!(df_healthy,history.loadhist["Healthy"],tend)

    CSV.write(folder_abs_path * "/LoadHist_Ill" * name_ext * ".csv", df_ill)
    CSV.write(folder_abs_path * "/LoadHist_Healthy" * name_ext * ".csv", df_healthy)

    return df_ill, df_healthy
end

function loadpostodataframe(history, folder_abs_path, name_ext = "")
    tend = history.par.historylength
    df_ill = empty_dataframe(history.par.Nloci,tend,zeros(Integer,2))
    df_healthy = empty_dataframe(history.par.Nloci,tend,zeros(Integer,2))

    populate_df_loadpos!(df_ill,history.loadpos["Ill"],tend)
    populate_df_loadpos!(df_healthy,history.loadpos["Healthy"],tend)

    CSV.write(folder_abs_path * "/LoadPos_Ill" * name_ext * ".csv", df_ill)
    CSV.write(folder_abs_path * "/LoadPos_Healthy" * name_ext * ".csv", df_healthy)

    return df_ill, df_healthy
end

empty_dataframe(ncols,nrows,zero::Number,startat=1) = DataFrame([Symbol(i) => [zero for _ in 1:nrows] for i in startat:startat+ncols-1])
empty_dataframe(ncols,nrows,zero::Vector,startat=1) = DataFrame([Symbol(i) => [copy(zero) for _ in 1:nrows] for i in startat:startat+ncols-1])
historylength(history) = length(history.mlp["PopSize"])

function populate_df_loadhist!(df,lh,tend)
    for t in 1:tend
        for (class, n) in zip(findnz(lh[t])...)
            df[t,class] += n
        end
    end
end

function populate_df_loadpos!(df,lh,tend) where {S<:Number,T<:Number}
    for t in 1:tend
        populate_df_loadpos_t!(df,lh[t],t)
    end
end

function populate_df_loadpos_t!(df,lp::Vector{SparseVector{S,T}},t) where {S<:Number,T<:Number}
    populate_df_loadpos_t!(df,lp[1],t,1)
    populate_df_loadpos_t!(df,lp[2],t,2)
end

function populate_df_loadpos_t!(df,lp::SparseVector,t,i)
    for (class,n) in zip(findnz(lp)...)
        df[t,class][i] += n
    end
end

function populate_df_loadpos_t!(df,lp::SparseVector,t)
    for (class,n) in zip(findnz(lp)...)
        df[t,class] += n
    end
end
