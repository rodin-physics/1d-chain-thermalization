filenames = readdir("data/Non_Thermal_Multi/")
filenames = filter(x -> first(x) !== '.', filenames)

root = joinpath(pwd(), "data/Non_Thermal_Multi/")

function check_str2(a)
    return tryparse(Float64, a) !== nothing
end


for name in filenames
    if isfile(joinpath(root, name)) && check_str2(string(first(name)))

        str = filter(x -> occursin("speed0",x), split(name,"_"))[1]
        speed = strip(str, [x for x in "speed0"])

        if isdir(joinpath(root, "Speed"*speed))
            mv(joinpath(root, name), joinpath(root, "Speed"*speed, name))
        else
            mkdir(joinpath(root, "Speed"*speed))
            mv(joinpath(root, name), joinpath(root, "Speed"*speed, name))
        end
    end
end
