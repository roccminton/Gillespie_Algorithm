using CairoMakie
using Random

include("/home/larocca/github/Gillespie_Algorithm/MainFunctions.jl")
import .Gillespie


f = Figure(resolution=(800,400))

N = 100

n_gen = 50
ε = 1/3    #percentage of space between generations compared to livespan

N0 = ceil(Integer,N/10)
t_grow = ceil(Integer,n_gen/4)

r = 1/500

PopSize(t) = t≥t_grow ? ceil(Integer,N/(1+((N-N0)/N0)*exp(-r*N*(t-t_grow)))) : N0

#---
#Wright Fisher
t = (1/n_gen) * (1-ε)
e = (1/n_gen) * ε

ax_wf = Axis(f[1,1],xlabel=L"Wright-Fisher")
hidedecorations!(ax_wf,label=false)
hidespines!(ax_wf)

for i in 1:n_gen
    vlines!(
        ax_wf,
        floor(Integer,(N-PopSize(i))/2):ceil(Integer,(N+PopSize(i))/2),
        ymin=(i-1)*(t+e),ymax=(i*t)+(i-1)*e,
        color=:blue
        )
end

#---
#Adaptive Dynamics
b = 1.0
d = 0.9
c(t) = t ≤ t_grow ? (b-d)/N0 : (b-d)/N

Nε = 5*ceil(Int,sqrt(N))

birthrate(N) = N*b
deathrate(N,t) = N*(d + N*c(t))
sort_free_index!(v) = sort!(v,by=x->abs(x-N/2-Nε),rev=true)

function birth!(PS,t)
    #Select free Index
    if !isempty(PS["Free_Index"])
        index = pop!(PS["Free_Index"])
        push!(PS["PopState"][index],t)
    else
        index = length(PS["PopState"])+1
        PS["N_max"] += 1
        push!(PS["PopState"],[t])
    end
    #Add Individual
    PS["PopSize"] += 1
    push!(PS["Occupied_Index"],index)
    #Reshuffle living Individuals
    shuffle!(PS["Occupied_Index"])
end

function death!(PS,t)
    #Select and individual at random
    index = pop!(PS["Occupied_Index"])
    #Delet Individal from population
    PS["PopSize"] -= 1
    push!(PS["PopState"][index],t)
    push!(PS["Free_Index"],index)
    #Resort free index
    sort_free_index!(PS["Free_Index"])
end

function execute!(PS,i,t)
    if i == 1
        birth!(PS,t)
    elseif i == 2
        death!(PS,t)
    else
        error("Unknown event with index $i.")
    end
end

function rate!(PS,t)
    PS["rates"][1] = birthrate(PS["PopSize"])
    PS["rates"][2] = deathrate(PS["PopSize"],t)
end

N₋ = ceil(Integer,N/2+Nε-N0/2)
N₊ = ceil(Integer,N/2+Nε+N0/2)

# PS = Dict(
#     "PopSize" => N0,
#     "PopState" => vcat(
#         [Float64[] for _ in 1:N₋],
#         [[0.0] for _ in N₋+1:N₊],
#         [Float64[] for _ in N₊+1:N+2Nε]
#         ),
#     "Free_Index" => sort_free_index!(vcat(collect(1:N₋),collect(N₊+1:N+2Nε))),
#     "Occupied_Index" => shuffle!(collect(N₋+1:N₊)),
#     "N_max" => N+2Nε,
#     "rates" => [birthrate(N0),deathrate(N0,0)]
# )
#
# t = 0
#
# while t <= n_gen
#     rate!(PS,t)
#     i, dt = Gillespie.nexteventandtime(PS["rates"])
#     iszero(i) && error("Population got extinct. Try again.")
#     #update time
#     global t += dt
#     execute!(PS,i,t)
# end
#
# #Kill all living individuals at end of simulation
# for i in PS["Occupied_Index"]
#     push!(PS["PopState"][i],n_gen)
# end

ax_ibm = Axis(
    f[1,3],
    xlabel=L"Adaptive Dynamics",
    yaxisposition=:right
    )
    hidedecorations!(ax_ibm,label=false)
    hidespines!(ax_ibm)


for x in 1:length(PS["PopState"])
    vlines!(
        fill(x,ceil(Integer,length(PS["PopState"][x])/2)),
        ymin=map(
            (x) -> x[2],
            Iterators.filter(
                ((x) -> isodd(x[1])),
                enumerate(PS["PopState"][x])
                )
            ) ./ n_gen,
        ymax=map(
            (x) -> x[2],
            Iterators.filter(
                ((x) -> iseven(x[1])),
                enumerate(PS["PopState"][x])
                )
            ) ./ n_gen,
        color=:orange
    )
end

#---
#Adjust ylim for Wright Fisher
n_max = PS["N_max"] - count(x->iszero(length(x)),PS["PopState"])

e = (n_max-N)/2

xlims!(ax_wf,low = -e, high = N+e)

#---
#Arrow for time
ax_arr = Axis(f[1,2],ylabel=L"Time")
hidedecorations!(ax_arr,label=false)
hidespines!(ax_arr)

up = 10
w = 1
h = 1


vlines!(ax_arr,0,ymin=0,ymax=up/(up+h),color=:black)
poly!(ax_arr,Point2f[(-w/2,up),(w/2,up),(0,h+up)],color=:black)
ylims!(ax_arr,low=0,high=up+h)

colsize!(f.layout,2,25)
colgap!(f.layout,1)

#---
#T-Growth line

# ax_dash = Axis(f[1,:])
# hidedecorations!(ax_dash)
# hidespines!(ax_dash)
#
# ylims!(ax_dash,high=n_gen)
# hlines!([t_grow],color=:red,linestyle=:dash)

#---

abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/"
save(abs_path * "CompareOverlappingvsNonOverlappingGensPopGrow_N=$N,tend=$n_gen.pdf", f, pt_per_unit=2)

f
