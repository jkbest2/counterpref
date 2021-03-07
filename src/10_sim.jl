using FisherySim
using SpatQSims
using Statistics
using StatsBase

import SpatQSims:
    survey_targeting,
    comm_catchability,
    sim_value,
    sim_values,
    simstudy_dir,
    simstudy_prefix


struct CounterPrefSpec{T} <: SpatQSimSpec
    realization::Int
    n_stations::Int
    init_pop::PopState{T}
    prep_file::String

    function CounterPrefSpec(realization::Int,
                             n_stations::Int,
                             init_pop::PopState{T},
                             prep_file::String = "prep.h5") where {T}
        new{T}(realization, n_stations, init_pop, prep_file)
    end
end

function CounterPrefSpec(realization::Int,
                         n_stations::Int,
                         prep_file::String = "prep.h5")
    prep = SpatQSimPrep(realization, prep_file)
    p0 = init_pop(prep)
    CounterPrefSpec(realization,
                    n_stations,
                    p0,
                    prep_file)
end


sim_value(spec::CounterPrefSpec) = spec.n_stations

# Choose survey stations with lower abundance
function survey_targeting(spec::CounterPrefSpec,
                          domain::AbstractFisheryDomain = domain(spec))
    survey_stations = vec(LinearIndices(domain.n)[3:5:98, 3:5:98])
    pop = spec.init_pop.P

    bin_pop = [sum(pop[i:(i + 4), j:(j + 4)]) for i in 1:5:96, j in 1:5:96]
    # Reverse the weighting
    inv_wt = 2mean(bin_pop) .- bin_pop
    # Make all weights nonnegative
    inv_wt .= inv_wt .- minimum(inv_wt)
    # Add a constant so that highest abundance stations have a non-zero
    # probability of selection
    inv_wt .+= 0.1 * maximum(inv_wt)

    stations = sample(survey_stations,
                      weights(inv_wt),
                      spec.n_stations;
                      replace = false,
                      ordered = true)

    FixedTargeting(stations)
end

# No spatially varying catchaiblity in this case
comm_catchbility(spec::CounterPrefSpec) = base_catchability()

sim_values(::Type{<:CounterPrefSpec}) = [50, 100, 200, 300, 400]

simstudy_dir(::Type{<:CounterPrefSpec}) = "sims/counterpref"
simstudy_prefix(::Type{<:CounterPrefSpec}) = "counterpref_"

repls = 1:25
prepf = "/home/jkbest/gscratch/spatq_sims/prep.h5"

run_sims(CounterPrefSpec, repls; prep_file = prepf, checkpoint = false)
