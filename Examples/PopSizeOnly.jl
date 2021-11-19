include("../Onetype.jl")

import .OneType
using Plots

#Setup Model parameters

t = 0:150

# logistic population growth
log_par = (
        birth = 1.0,
        death = 0.9,
        competition = 10^(-5)
        )
logistic = OneType.modelsetup("Logistic")
loghis = OneType.rungillespie(t,100,log_par,logistic)

# logistic population growth rescaled
log_par_res = (log_par...,K = 100)
logistic_res = OneType.modelsetup("Logistic",rescaled=true)
logreshis = OneType.rungillespie(t,100.0,log_par_res,logistic_res)

#exponential growth
yule_par = (birth=0.05,)
yule = OneType.modelsetup("Yule")
yulehis = OneType.rungillespie(t,10,yule_par,yule)

#linear birht and death process
linear_par = (birth=0.5,death=0.5)
linear = OneType.modelsetup("Linear")
linearhis = OneType.rungillespie(t,5000,linear_par,linear)

#birth death process with immigration
imm_par = (linear_par...,immigration = 1)
imm = OneType.modelsetup("Immigration")
immhis = OneType.rungillespie(t,5000,imm_par,imm)

#----Plotting---

tfs = 10

plog = plot(
        t,[loghis,logreshis],
        label=["Logistic growth" "Logistic growth rescaled"],
        legend=:bottomright,
        title = "Logistic Population Growth",
        titlefontsize = tfs
)

pyule = plot(
        t,yulehis,
        title="Exponential Growth",
        color=:red, legend=false,
        titlefontsize = tfs
        )
plinear = plot(
        t,linearhis,
        title="Linear Birth and Death Process",
        color=:darkgreen,
        legend=false,
        titlefontsize = tfs
        )
pimm = plot(
        t,immhis,
        title = "Linear Process with Immigration",
        color=:darkorange,
        legend=false,
        titlefontsize = tfs
        )

plot(plog,pyule,plinear,pimm, layout=(2,2),size=(800,600))
