using Random
include("../src/main.jl")
## System Parameters
d = 60
ωmax = 10
τmax = 200
n_masses = 1000

δ = (1 / ωmax) / d              # Time step
n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
τs = δ * (1:n_pts)

## Plotting homogeneous motion of single chain atom
# Ignore hidden files
ignorefile(str) = (str[1] == '.')

root = "precomputed/rH"
# filenames = joinpath.(root, filter(!ignorefile, readdir(root)))

filenames = ["precomputed/rH/rH_ωmax10_d60_ωT1.0e-5_τ200.jld2",
            "precomputed/rH/rH_ωmax10_d60_ωT1.0_τ200.jld2",
            "precomputed/rH/rH_ωmax10_d60_ωT10.0_τ200.jld2",
            "precomputed/rH/rH_ωmax10_d60_ωT100.0_τ200.jld2",
            "precomputed/rH/rH_ωmax10_d60_ωT500.0_τ200.jld2"]

colors = reverse([my_red, my_vermillion, my_green, my_blue, my_black])

idx = 1    # Index of chain atom

data = Array{Float64}[]
ωTs = Float64[]
for ind in 1:length(filenames)
    file = load_object(filenames[ind])
    ρs = file.ρHs
    push!(data, ρs[idx , :])
    push!(ωTs, file.ωT)
end

perm = sortperm(ωTs)
data = data[perm]

fig = Figure(resolution = (1800, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 50)
ax = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
for ind in reverse(1:length(data))
    σ = data[ind]
    lines!(ax, τs, σ, linewidth = 4, color = colors[ind], label = L"\omega_T = %$(ωTs[ind])")
end

axislegend(ax)
CairoMakie.xlims!(ax, 0, 200)
fig

## Check standard deviation

qs = range(0, π / 2, length = round(n_masses/2) |> Integer)
ωs = ω.(ωmax, qs)

Random.seed!(150)
ϕs = 2 * π * rand(length(qs))
ωTs = range(1e-5, 2500, length = 300)

# Analytic standard deviation of chain atom homogeneous trajectory
function analytic_std_dev(ωT)
    int_func(x) = coth(ω(ωmax, x) / (2 * ωT)) / ω(ωmax, x)
    return sqrt((1/π) * quadgk(int_func, 0, π/2)[1])
end
std_devs = map(analytic_std_dev, ωTs)

# Zero-point motion
std_dev_T0 = sqrt((1/(π * ωmax)) * ellipk(1 - (1 / ωmax^2)))

# Calculate standard deviation of homogeneous trajectory
dev = Float64[]
@showprogress for ωT in ωTs
    ζs = ζq.(ωs, ωT)
    ind = 0

    res = Statistics.stdm(map(n -> ζH(n, δ, ζs, ϕs, ωs, qs, ind), 1:n_pts), 0)

    push!(dev, res)
end

fig2 = Figure(resolution = (1600, 1200), font = "CMU Serif", fontsize = 36, figure_padding = 50)
ax1 = Axis(fig2[1, 1], xlabel = L"\omega_T", ylabel = "Std Dev")

lines!(ax1, ωTs, std_devs, color = my_blue, linewidth = 4)
hlines!(ax1, [std_dev_T0], color = my_black, linewidth = 4, linestyle = :dashdot)
scatter!(ωTs, dev, marker = :circle, color = my_sky, markersize = 20)

CairoMakie.xlims!(ax1, 0.0, 2500.0)
CairoMakie.ylims!(ax1, 0.0, nothing)
fig2
