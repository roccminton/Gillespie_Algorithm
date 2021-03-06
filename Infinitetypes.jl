#=
    module InfiniteTypes

    Implementation of a Logistic Model with a infinite number of possible types.
    Mutation are possible at an independent rate μ and are specified by a
    transition kernel from one type to the next.

    Rates are represended as a Dictionary of Vectors, where the frist entry
    of the Vector gives the birh and the second entry the death rate of the
    individual type with that key.
=#

module InfiniteTypes

include("MainFunctions.jl")
import .Gillespie

using NamedTupleTools #for fieldnames

export rates!, execute!

function rungillespie(time,x0,n_x0,par;rescaled=false)
    ps0 = initialpopulation(x0,n_x0,par)
    his = setuphistory(time,ps0)

    Gillespie.run_gillespie!(
        time,ps0,
        setupparameter(ps0,par,rescaled),
        execute!,rates!,
        setuprates(par.birth,par.death,ps0),
        his
        )

    return his
end

function generaterescaledparam(par)
    #generate rescaled functions
    birth = x -> par.K * par.birth(x)
    death = x -> par.K * par.death(x)
    competition = (x,y) -> par.K * par.competition(x,y)

    #copy other fields
    fields = (field for field in fieldnames(par) if field ∉ [:birth,:death,:competition])
    values = (par[field] for field in fields)

    #return rescaled tuple with rescaled mutation and rescaled diff
    return NamedTuple{(fields...,:diff,:birth,:death,:competition)}(
        values...,
        inv(par.K),
        birth,
        death,
        competition)
end

function setupparameter(ps0,par,rescaled)
    with_comp_par = (
        par...,
        compdict = generatecompdict(ps0,par.competition)
        )
    if rescaled
        return generaterescaledparam(with_comp_par)
    else
        return (with_comp_par...,diff=one(eltype(valtype(ps0))))
    end
end

function setuprates(b,d,x0)
    Individuals = collect(keys(x0))
    return Dict{
            keytype(x0),
            Vector{promote_type(typeof(b(Individuals[1])),typeof(d(Individuals[1])))}
        }()
end

initialpopulation(x,nₓ,par) = Dict(x => [nₓ,par.birth(x),par.death(x)])

function setuphistory(time,n₀)
    T = eltype(valtype(n₀))
    l = length(time)
    return Dict(x=>zeros(T,l) for x in keys(n₀))
end

function generatecompdict(ps,competition)
    IndividualType = keytype(ps)
    Individuals = collect(keys(ps))
    #generate empty dictionary
    C = Dict{
        Tuple{IndividualType,IndividualType},
        typeof(competition(Individuals[1],Individuals[1]))
        }()
    #populate dictionary
    for x in keys(ps), y in keys(ps)
        C[(x,y)] = competition(x,y)
    end
    return C
end

function rates!(rates::Dict,ps::Dict,pr)
    for (x,vₓ) in ps
        nₓ = vₓ[1]
        !haskey(rates,x) && (rates[x] = valtype(rates)(undef,2))
        #birthrate n_x * b(x)
        rates[x][1] = nₓ*vₓ[2]
        #deathrate n_x * (d(x) + Σ c(x,y) n_y)
        rates[x][2] = nₓ* vₓ[3]
        for (traittuple,c) in pr.compdict
            t₁,t₂ = traittuple
            t₁ == x && (rates[x][2] += nₓ * ps[t₂][1] * c)
        end
    end
end

function birth!(ps, rates, pr, trait)
    #Birth with or without mutation
    if pr.μ > 0.0 && rand() ≤ pr.μ
        #mutate to new type/species and add to species
        new_trait = pr.mutate(trait)
        #setup the size of the new type
        if haskey(ps,new_trait)
            ps[new_trait][1] += pr.diff
        else
            addnewtrait!(ps,rates,pr,new_trait)
        end
    else
        ps[trait][1] += pr.diff
    end
    nothing
end

function addnewtrait!(ps,rates,pr,trait)
    #add to population state
    ps[trait] = [pr.diff,pr.birth(trait),pr.death(trait)]
    #set competition
    for other_trait in keys(ps)
        pr.compdict[(trait,other_trait)] = pr.competition(trait,other_trait)
        pr.compdict[(other_trait,trait)] = pr.competition(other_trait,trait)
    end
end

function death!(ps,trait,pr)
    ps[trait][1] -= pr.diff
end

function execute!(i,trait,ps,rates,pr)
    if i==1
        birth!(ps, rates, pr, trait)
    elseif i==2
        death!(ps,trait,pr)
    else
        error("Index Error: Unknown event #$i")
    end
end

monoeq(x,par) = (par.birth(x) - par.death(x)) / par.competition(x,x)

end  # end of module InfiniteTypes
