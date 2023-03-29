include("../../src/main.jl")

## Parameters
α = 40
μ = 1
ωmax = 10
nChain = 3
Φ0 = 2.0
λ = 4.0
# Load memory kernel and thermal trajectory
system = load_object("precomputed/systems/System_ωmax10_d60_l10000_τmax5.jld2")

## Computation 
σdot0 = 50.0
ωT = 0.0
colors = [my_blue, my_green, my_vermillion]
ωTs = [0.0, 5.0, 25.0]




## Plotting density compared to predicted distribution
fig = Figure(resolution = (1200, 1000), fontsize = 40, figure_padding = 30)
ax1 = Axis(fig[1, 1], xlabel = L"\Delta", ylabel = "Probability Density", xgridvisible = false, ygridvisible = false)

for ωT in ωTs 
    Δs = readdlm("data/Thermal/ωT$(ωT)/delta_ωmax10_speed$(σdot0)_ωT$(ωT).dat") |> vec
    hist!(ax1, Δs, npoints = 200, label = "Numerical", color = (colors[findfirst(x -> x == ωT, ωTs)],  0.5), normalization = :pdf, bins = 100)

    pred_mean = Δ_thermal_analytic(σdot0, Φ0, λ, ωmax, ωT)
    pred_var = Δ_thermal_variance(σdot0, Φ0, λ, ωmax, ωT)

    density!(ax1, rand(Normal(pred_mean, √(pred_var)), 1000000), color = :transparent, strokearound = true, strokewidth = 4, label = "Predicted", strokecolor = colors[findfirst(x -> x == ωT, ωTs)])

end

temp_color = [PolyElement(color = (color, 0.5), strokecolor = :transparent) for color in colors]
res_type = [PolyElement(color = (my_black, 0.5), strokecolor = :transparent), LineElement(color = (my_black, 0.5), linewidth = 4)]


Legend(fig[1,1],
    [temp_color, res_type],
    [string.(ωTs), ["Numerical", "Predicted"]],
    [L"\omega_T", "Results"],
    tellheight = false,
    tellwidth = false,
    halign = :right,
    valign = :top,
    margin = (10, 10, 10, 10)) 

# Legend(fig[1, 1],
#     [elem_1, elem_2, elem_3, elem_4, elem_5],
#     ["Line & Marker", "Poly & Line", "Line", "Marker", "Poly"],
#     tellheight = false,
#     tellwidth = false,
#     halign = :right,
#     valign = :top,
#     margin = (10, 10, 10, 10))


# axislegend(ax1, labelsize = 35, position = :rt)
xlims!(ax1, -2.5, 2.5)
ylims!(ax1, 0, nothing)
fig