
include("/home/larocca/github/Gillespie_Algorithm/ToolBoxPlotting.jl")


tend =15_000
K = 10_000
N = 90
dni = 0.4
i = 1

#Save simulations
#---
abs_path = "/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=$K,dni=$dni,N=$N"

#Load Data
#---

d0 = load(abs_path * "/truerates$i.jld")
d1 = load("/home/larocca/github/Gillespie_Algorithm/DiploidModel/Data/NoRecombination/K=10000,dni=0.4,N=90/r=1.jld")
d0["ill"] = d0["Ill"]
d1["ill"] = d1["Ill"]
#Plot Extinction of Truerates Simulations
#---
ext = findextime(d)
max_load = max(maximum(ml(d0)),maximum(ml(d1)))
max_prev = max(maximum(prev(d0)),maximum(prev(d1)))

f = Figure(size=(800,450),fontsize=17)

axload, axprev, loadline, prevline = add_mlp_plot!(
    f,d0,
    tend,
    c2 = RGBf(0/255,80/255,162/255),
    c1 = RGBf(198/255,38/255,6/255)
    #ceil(Integer,ext*1.1)
    )

ylims!(axload,(0,max_load))
ylims!(axprev,(0,max_prev))

#text!(ceil(Integer,tend*0.1),ceil(Integer,max_load*0.85),text=L"$r=0",fontsize=30,align=(:left,:center))
text!(ceil(Integer,tend*0.1),ceil(Integer,max_load*0.85),text=L"μ=%$(dni/2)",fontsize=30,align=(:left,:center))
text!(ceil(Integer,tend*0.1),ceil(Integer,max_load*0.75),text=L"N=%$(N)",fontsize=30,align=(:left,:center))

#save("/home/larocca/github/Diss/Defense/truerates1.pdf",f, pt_per_unit=2)

f
#---

#CovMatrix Gif Frame in Frame
#---
# cormatrixs = CovToCor.(d["CovMatrixList"])
#
# #generations per frame
# rate = 100
#
# t = Observable(1)
#
# f = Figure(resolution=4 .* (682,381))
#
# axload, axprev, loadline, prevline = add_mlp_plot!(
#     f,d,
#     tend,
#     1
#     )
#
# vlines!(axload,@lift([$t]),color=:black,linewidth=1)
#
# axcor = inset_axis!(f,axload,extent=(0.5,1,0,0.6),alignmode=Outside(10))
# heatmap!(axcor,@lift(cormatrixs[$t]),clim=(0,1))
#
# hidedecorations!(axcor)
#
# Colorbar(f[2,:], limits = (0, 1),vertical=false,tellheight=true,tickalign=1)
#
# rowgap!(f.layout,1)
#
# timestamps = 1:rate:tend+1
#
# record(f,abs_path * "/corrmatrix_rate=$rate.gif", timestamps) do x
#     t[]=x
# end

#CorrCluster Plot
#---

# f = Figure(resolution=1.75 .* (682,381))
#
# T = round.(Integer,0:10:tend)
#
# axload, axprev, loadline, prevline = add_mlp_plot!(
#     f,d,
#     tend,
#     1
#     )
#
# axcount = Axis(f[2,1],yaxisposition =:right,xticksvisible = false,yticklabelcolor = :green,ylabel="Genes in Cluster")
# axn = Axis(f[2,1],yticklabelcolor = :blue,ylabel="Number of Clusters",xlabel="Time")
#
# hidexdecorations!(axcount)
#
# lines!(axcount,T,cluster_count.(d["CorrCluster"][T.+1]),color=:green)
# lines!(axn,T,n_cluster.(d["CorrCluster"][T.+1]),color=:blue)
#
# xlims!(axcount,(0,tend))
# xlims!(axn,(0,tend))
# ylims!(axcount,low=0)
# ylims!(axn,low=0)
#
# axn.xticks = (
#     range(0, tend; length = 5),
#     [string.(round.(Integer, x / 1000)) .* "K" for x in range(0,tend;length=5)]
#      )
#
# #save(abs_path *"/fig_CorrCluster.pdf",f, pt_per_unit=2)
#
# f
