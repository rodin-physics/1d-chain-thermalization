include("../src/main.jl")
using Peaks
using DelimitedFiles

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), "data/final_tau800/")))

for ii in 1:length(filenames)
    println(ii, " out of ", length(filenames))
    data = load_object(joinpath(pwd(), "data/final_tau800/", filenames[ii]))
    (τs, speeds) = particle_speed(data)

    if data.Φ > 0 && all(speeds.>=0)
        pks, vals = findmaxima(speeds)
        cd("data/proc_tau800/")
        writedlm(filenames[ii] * ".dat", [τs[pks] vals])
        cd("../../")

    elseif data.Φ < 0 && all(speeds.>=0)
        pks, vals = findminima(speeds)
        cd("data/proc_tau800/")
        writedlm(filenames[ii] * ".dat", [τs[pks] vals])
        cd("../../")
    end
end
