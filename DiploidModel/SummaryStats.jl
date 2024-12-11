#Define Type for Population History
struct MLPHistory
    mlp :: Dict
    par :: NamedTuple
end

#Setup an empty population history to fill during the runtime of the simulation
function DiploidModel2.setup_pop_hist(par,n₀,l)
    history = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))
    return MLPHistory(history,par)
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
    #safe all the parameters as seperate entries
    for (k,v) in zip(keys(h.par),h.par)
        key = String(k)
        key ∈ safe_parameter && (safe_h[key] = v)
    end
    return safe_h
end

#overwrite the basic choice of the default saveonestep function in Gillespie if necessary
DiploidModel2.choosestatsfunction(population_history::MLPHistory) = saveonestep!

#and define your own saveonestep! function
saveonestep!(ph::MLPHistory,index,ps,par) = DiploidModel2.Gillespie.saveonestep!(ph.mlp,index,ps,par)

DiploidModel2.updatestats_death!(ps,par,index) = nothing
DiploidModel2.updatestats_birth!(ps,par,index) = nothing
