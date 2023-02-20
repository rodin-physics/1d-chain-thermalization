include("../../src/main.jl")

step_size = 20
## General Example

fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma")

# REPULSIVE

data = load_object(
    "data/Thermal/Single/Single_sigma0[180]_sigmadot0[120]_MemInf_lambda4_Phi20_mu1_d60_bias0.0_omegaT25.0_tau200.jld2",
)
δ = data.τs[2] - data.τs[1]
rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
# mx = maximum(abs.(rr))
mx = 1.0
hm = heatmap!(
    ax1,
    data.τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)
lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)
xlims!(ax1, (0, 200))
ylims!(ax1, (0, 10000))

# lines!(
#     ax1,
#     0:0.1:6,
#     π * 10 * 9 * (0:0.1:6) .+ data.σs[1][1],
#     color = my_black,
#     linewidth = 4,
#     linestyle = :dash,
# )
Colorbar(fig[1, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
# ATTRACTIVE

data = load_object("data/Thermal/Single/Single_sigma0[180]_sigmadot0[120]_MemInf_lambda4_Phi-20_mu1_d60_bias0.0_omegaT25.0_tau200.jld2")
δ = data.τs[2] - data.τs[1]
rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
# mx = maximum(abs.(rr))
mx = 1.0
hm = heatmap!(
    ax2,
    data.τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)
lines!(ax2, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)
xlims!(ax2, (0, 200))
# ylims!(ax2, (0, 10000))

# lines!(
#     ax2,
#     0:0.1:6,
#     π * 10 * 9 * (0:0.1:6) .+ data.σs[1][1],
#     color = my_black,
#     linewidth = 4,
#     linestyle = :dash,
# )
Colorbar(fig[2, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
fig