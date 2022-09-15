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

dir_name = "data/final_tau800/"
final_dir = "data/proc_tau800"

## Process data in dir_name and extract velocities between chain masses
filenames = filter(x -> first(x) !== '.', readdir(joinpath(pwd(), dir_name)))

for ii in 1:length(filenames)
    println(ii, " out of ", length(filenames))
    data = load_object(joinpath(pwd(), dir_name, filenames[ii]))
    (τs, speeds) = particle_speed(data)

    if data.Φ > 0 && all(speeds.>=0)
        pks, vals = findmaxima(speeds)
        cd(final_dir)
        writedlm(filenames[ii] * ".dat", [τs[pks] vals])
        cd("../../")

    elseif data.Φ < 0 && all(speeds.>=0)
        pks, vals = findminima(speeds)
        cd(final_dir)
        writedlm(filenames[ii] * ".dat", [τs[pks] vals])
        cd("../../")
    end
end
