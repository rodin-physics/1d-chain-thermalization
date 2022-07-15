include("../src/main.jl")

step_size = 10

# Speed of particle over time
function particle_speed(data)
    σs = data.σs |> vec
    τs = data.τs
    δ = τs[2] - τs[1]
    speeds = [((σs[ii+1] - σs[ii]) / δ) for ii in 1:(length(σs)-1)]

    return (τs[2:end], speeds)
end

## Plotting
fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma")

# REPULSIVE

data = load_object(
    "data/non_thermal/Single_σ0[200]_σdot0[50]_MemInf_λ4_Φ8_μ1_d300_bias0.0_ΩTnothing_τ100.jld2")

# Find all times where resonances occur
resonance_speed(n) = (2 * data.ωmax * data.α) / ((2 * n) + 1)
(all_τs, σdots) = particle_speed(data)
res_times = [findall(x -> isapprox(x, ii, atol = 0.01) == true, σdots) for ii in resonance_speed.(10:15)]

δ = data.τs[2] - data.τs[1]

τ_max = 100
idx = findall(data.τs .< τ_max)
τs = data.τs[idx]
rr = reduce(hcat, [data.ρs[ii, idx] .- ii * data.α for ii = 1:size(data.ρs)[1]])

mx = maximum(abs.(rr)) / 1
mx = 0.1
hm = heatmap!(
    ax1,
    τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)
lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)

# for ii in last.(res_times)
#     lines!(ax1, (0:0.1:30) .+ all_τs[ii], π * data.α * (data.ωmax - 1) * (0:0.1:30) .+ data.σs[ii], color = my_sky, linestyle = :dash, linewidth = 4)
#
#     lines!(ax1, (0:0.1:10) .+ all_τs[ii], -π * data.α * (data.ωmax - 1) * (0:0.1:10) .+ data.σs[ii], color = my_sky, linestyle = :dash, linewidth = 4)
# end

lines!(
    ax1,
    0:0.1:10,
    π * data.α * (data.ωmax - 1) * (0:0.1:10) .+ data.σs[1][1],
    color = my_black,
    linewidth = 4,
    linestyle = :dash,
)
Colorbar(fig[1, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
xlims!(ax1, (0, 100))
ylims!(ax1, (0, 4800))
# # ATTRACTIVE

data = load_object(
    "data/non_thermal/Single_σ0[220]_σdot0[50]_MemInf_λ4_Φ-8_μ1_d300_bias0.0_ΩTnothing_τ100.jld2",
)
# Find all times where resonances occur
resonance_speed(n) = (2 * data.ωmax * data.α) / ((2 * n) + 1)
(all_τs, σdots) = particle_speed(data)
res_times = [findall(x -> isapprox(x, ii, atol = 0.004) == true, σdots) for ii in resonance_speed.(8:15)]

δ = data.τs[2] - data.τs[1]
τ_max = 100
idx = findall(data.τs .< τ_max)
τs = data.τs[idx]
rr = reduce(hcat, [data.ρs[ii, idx] .- ii * data.α for ii = 1:size(data.ρs)[1]])
mx = 0.1
# mx = maximum(abs.(rr)) / 3
hm = heatmap!(
    ax2,
    τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)
lines!(ax2, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)

# for ii in first.(res_times)
#     lines!(ax2, (0:0.1:30) .+ all_τs[ii], π * data.α * (data.ωmax - 1) * (0:0.1:30) .+ data.σs[ii], color = my_sky, linestyle = :dash, linewidth = 4)
#
#     lines!(ax2, (0:0.1:10) .+ all_τs[ii], -π * data.α * (data.ωmax - 1) * (0:0.1:10) .+ data.σs[ii], color = my_sky, linestyle = :dash, linewidth = 4)
# end

lines!(
    ax2,
    0:0.1:10,
    π * data.α * (data.ωmax - 1) * (0:0.1:10) .+ data.σs[1][1],
    color = my_black,
    linewidth = 4,
    linestyle = :dash,
)
Colorbar(fig[2, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)

xlims!(ax2, (0, 100))
ylims!(ax2, (0, 4800))
fig
# save("General_Example.pdf", fig)
