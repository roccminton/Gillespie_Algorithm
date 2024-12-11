include("DiploidModel2.jl")
import .DiploidModel2

#choose birht death functions package
include("/home/larocca/github/Gillespie_Algorithm/DiploidModel/BirthDeathwithRec.jl")
#choose stats function package
include("/home/larocca/github/Gillespie_Algorithm/DiploidModel/AllStatswithHaploidMutClasses.jl")


#---
using JLD
using CairoMakie

#---
healthy_pop(K) = Dict(
    "PopSize" => K,
    "Ill" => 0,
    "ML" => 0
    )

uniform_competition(b,d,K) = (b-d)/K

model_parameter(K,dni,N,b,d,c,rec,birthrates,path) = (
        birth = b,
        death = d,
        competition = c,
        μ = dni,
        Nloci = N,
        K = K,
        rates = birthrates,
        recombination = rec,
        path = path
)

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)
prev(h) = replace_NaN(h.mlp["Ill"] ./ h.mlp["PopSize"])
ml(h) = replace_NaN(h.mlp["ML"] ./ h.mlp["PopSize"])

#---
function dosimulation(N,dni,filename;K=10_000,tend=100_000,b=1.0,d=0.9,c=nothing,rec=0,birthrates="allbirthrates!",nruns=3)
    abs_path = "/media/larocca/PortableSSD/Data/NoRecombination"
    abs_path = mkpath(abs_path)
    filename = "/" * filename

    isnothing(c) && (c=uniform_competition(b,d,K))

    println("Currently at dni = $(round(dni,digits=2)), N = $N")
    f = Figure(;size=(600,300*nruns))
    for i in 1:nruns
        #Run Simulation
        h = DiploidModel2.rungillespie(
            1:tend+1,
            healthy_pop(K),
            model_parameter(K,dni,N,b,d,c,rec,birthrates, abs_path*filename*"_$i")
        )
        #Save Data
        save(abs_path * filename *"_$i.jld",convertforsaving(h))
        #Create Plot
        add_mlp_plot!(f,h,i,tend)
        #free space after safing the file
        h = nothing
    end
    #save data and plot
    save(abs_path * filename *"_overview.pdf", f)
end

function add_mlp_plot!(f,h,i,tend)
    #create two axis which share a x-axis
    axprev = Axis(
        f[i, :],
        yticklabelcolor = :orange,
        yaxisposition = :right,
        xaxisposition = :top,
        xticksvisible = false,
        ytickformat = x -> string.(round.(Integer, x * 100)) .* "%",
        ylabel = "Prevalence",
    )

    axload = Axis(
        f[i, :],
        yticklabelcolor = :red,
        ylabel = "Mutation Burden",
        xlabel = "Time",
    )

    axload.xticks = (
        range(0, tend; length = 5),
        [
            string.(round.(Integer, x / 1000)) .* "K" for
            x in range(0, tend; length = 5)
        ],
    )

    hidespines!(axprev)
    hidexdecorations!(axprev)

    #set the xlims
    linkxaxes!(axload, axprev)
    xlims!(axload, (0, tend))
    xlims!(axprev, (0, tend))

    #plot the data
    prevline = lines!(
        axprev,
        0:tend,
        prev(h)[1:tend+1],
        color = :orange,
        label = "Prevalence",
    )
    loadline =
        lines!(axload, 0:tend, ml(h)[1:tend+1], color = :red, label = "Mutation Burden")
end

#---
