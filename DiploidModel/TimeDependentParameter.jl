

#General Model Functions for parameter changes during runtime

function changeparameter(par,newpar)
    return merge(par,newpar)
end

function runseveralgillespie(totaltime,switch::Vector,parameter_changes,n₀,model_parameter)

    #setup empty population history
    l = length(totaltime)
    population_history = Dict(x=>zeros(valtype(n₀),l) for x in keys(n₀))

    #setup empty rates vector
    initrates = Vector{typeof(model_parameter.birth)}(undef,2)

    #choose rates function
    if haskey(model_parameter,:rates)
        rates! = getfield(DiploidModel2,Symbol(model_parameter.rates))
    else
        rates! = truerates!
    end

    par = setupparameter(model_parameter,n₀,l)

    #save initial population
    for (x,nₓ) in n₀
        #add the population history for the given trait
        population_history[x][1] = nₓ
    end

    #switching times
    times = timesfromswitches(totaltime,switch)

    #execute simulation
    for (i,t) in enumerate(times)
        Gillespie.run_gillespie!(
            t,n₀,
            par,
            execute!,
            rates!,
            initrates,
            population_history,
            hstart = t[1]
            )
        #change parameter for next block, unless no block follows
        i < length(times) && (par = changeparameter(par,parameter_changes[i]))
    end
    return population_history
end

function runseveralgillespie(totaltime,switch::Number,par_changes,x0,parameters)
    return runseveralgillespie(totaltime,[switch],par_changes,x0,parameters)
end
