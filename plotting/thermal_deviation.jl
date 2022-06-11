include("../src/main.jl")

# data_rep = load_object("data/Thermal/Single_σ0[35]_σdot0[20]_MemInf_λ1_Φ1_μ1_d60_ΩT1.0e-6_τ120.jld2")
#
# data_attr = load_object("data/Thermal/Single_σ0[35]_σdot0[20]_MemInf_λ1_Φ-1_μ1_d60_ΩT1.0e-6_τ120.jld2")
#
# α = data_rep.α
# τs = data_rep.τs
# δ = τs[2] - τs[1]
#
# ρs_rep = data_rep.ρs
# ρs_attr = data_attr.ρs
#
# idx = 150    # Index of chain atom
#
# fig = Figure(resolution = (1200, 1800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
# ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma")
#
# lines!(ax1, τs, ρs_rep[idx, :] .- (idx * α), linewidth = 5)
# lines!(ax2, τs, ρs_attr[idx, :] .- (idx * α), linewidth = 5)
#
# CairoMakie.xlims!(ax1, 0, 10)
# CairoMakie.xlims!(ax2, 0, 10)
# fig

## Check standard deviation
ωmax = 10
n_masses = 200              # Number of chain masses for simulating ρ0
qs = range(0, π / 2, length = round(n_masses / 2) |> Integer)
ωs = ω.(ωmax, qs)
Random.seed!(150)
ϕ1s = 2 * π * rand(length(qs))
Random.seed!(152)
ϕ2s = 2 * π * rand(length(qs))
ωTs = range(1e-5, 10, length = 30)

# Analytic standard deviation of chain atom homogeneous trajectory
function analytic_std_dev(ωT)
    int_func(x) = coth(ω(ωmax, x) / (2 * ωT)) * (1 / ω(ωmax, x))
    return (1/π) * quadgk(int_func, 0, π/2)[1]
end

std_dev_T0 = (1/(π * ωmax)) * ellipk(1 - (1 / ωmax^2))
std_devs = map(analytic_std_dev, ωTs)

dev = Float64[]
for ωT in ωTs
    println("ωT is ", ωT)
    ζs = ζq.(ωs, ωT)
    idx = 10
    res = Statistics.stdm(map(n -> (ζH_sin(n, δ, ζs, ϕ1s, ωs, qs, idx) + ζH_cos(n, δ, ζs, ϕ2s, ωs, qs, idx)), 1:n_pts), idx)
    push!(dev, res .- idx)
end

fig = Figure(resolution = (1600, 1600), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1, 1], xlabel = L"\omega_T", ylabel = "Std Dev")
scatter!(ωTs, dev, marker = :circle, markersize = 20)
lines!(ax1, ωTs, std_devs, color = :blue, linewidth = 4)
hlines!(ax1, [std_dev_T0], color = :black, linewidth = 4, linestyle = :dashdot)
fig
