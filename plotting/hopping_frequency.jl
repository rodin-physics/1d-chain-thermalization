include("../src/main.jl")

root = "Data/Thermal"
colors = reverse([my_red, my_vermillion, my_green, my_sky, my_blue, my_black])
## Load zero-speed single-particle data
rep_files = [
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ1_μ1_d60_ΩT1.0e-5_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ1_μ1_d60_ΩT1.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ1_μ1_d60_ΩT50.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ1_μ1_d60_ΩT100.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ1_μ1_d60_ΩT250.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ1_μ1_d60_ΩT500.0_τ100.jld2",
]

attr_files = [
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ-1_μ1_d60_ΩT1.0e-5_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ-1_μ1_d60_ΩT1.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ-1_μ1_d60_ΩT50.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ-1_μ1_d60_ΩT100.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ-1_μ1_d60_ΩT250.0_τ100.jld2",
    "Single_σ0[225]_σdot0[0]_MemInf_λ1_Φ-1_μ1_d60_ΩT500.0_τ100.jld2",
]

rep_data = load_object.(joinpath.(root, rep_files))
attr_data = load_object.(joinpath.(root, attr_files))

# Mean squared displacement of particle
function mean_squared_disp(σs::Vector{Float64})
    vec1 = σs[2:end]
    vec2 = σs[1:end-1]
    return msd(vec1, vec2)
end


# Determine avg time taken for mobile particle to hop to next chain atom pair
function hopping_freqs(sol::SystemSolution)
    σs = sol.σs
    τs = sol.τs
    ρs = sol.ρs
    σ0 = σs[1]
    # Initialise array
    hops = Int64[]

    for ind = 1:(length(τs)-1)
        # Determine position of mobile particle in relation to chain
        curr = searchsortedfirst(ρs[:, ind], σs[ind])
        nxt = searchsortedfirst(ρs[:, (ind+1)], σs[ind+1])

        # Add the difference between initial and final positions
        push!(hops, abs(nxt - curr))
    end

    # Find the time between movements
    indices = findall(x -> x != 0, hops)
    times = [(τs[indices[ii+1]] - τs[indices[ii]]) for ii = 1:(length(indices)-1)]

    # Return the time between movements and when those movements occured
    return (τs[indices[2:end]], times)
end


## Plotting
fig = Figure(
    resolution = (1200, 1600),
    font = "CMU Serif",
    fontsize = 36,
    figure_padding = 30,
)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = "Time between Hops")
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = "Time between Hops")

for ind = 1:length(rep_files)
    (xs, ys) = hopping_freqs(rep_data[ind])
    scatter!(
        ax1,
        xs,
        ys,
        color = colors[ind],
        markersize = 8,
        label = L"\omega_T = %$(rep_data[ind].ωT)",
    )
    hlines!(ax1, [mean(ys)], color = colors[ind], linewidth = 5)
end

for ind = 1:length(attr_files)
    (xs, ys) = hopping_freqs(attr_data[ind])
    scatter!(
        ax2,
        xs,
        ys,
        color = colors[ind],
        markersize = 8,
        label = L"\omega_T = %$(attr_data[ind].ωT)",
    )
    hlines!(ax2, [mean(ys)], color = colors[ind], linewidth = 5)
end
axislegend(ax1)
CairoMakie.xlims!(ax1, 0.0, 100)
CairoMakie.ylims!(ax1, 0.0, 1.2)

axislegend(ax2)
CairoMakie.xlims!(ax2, 0.0, 100)
CairoMakie.ylims!(ax2, 0.0, 3.2)
fig
