
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

historytodataframe(history) = DataFrame(
    PopSize=history["PopSize"],
    Ill=history["Ill"],
    Mutation=history["ML"],
    Ccuts = vcat(history["cutsat"],fill(NaN,length(history["PopSize"])-length(history["cutsat"])))
    )

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

function maxmutatioins(hhist,key="both")
    if key=="both"
        return max(findnz(hhist["Ill"][end])[1][end],findnz(hhist["Healthy"][end])[1][end])
    else
        return findnz(hhist[key][end])[1][end]
    end
end
