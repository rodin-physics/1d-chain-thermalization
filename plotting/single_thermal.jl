include("../src/main.jl")

step_size = 20
## General Example

fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 36)
ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma")

# REPULSIVE

data = load_object(
    "data/Thermal/Single_σ0[35]_σdot0[20]_MemInf_λ1_Φ1_μ1_d60_ΩT1.0e-6_τ120.jld2",
)
δ = data.τs[2] - data.τs[1]
rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
mx = maximum(abs.(rr))
hm = heatmap!(
    ax1,
    data.τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)
lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)
# xlims!(ax1, (0, 80))
ylims!(ax1, (0, 2000))

lines!(
    ax1,
    0:0.1:6,
    π * 10 * 9 * (0:0.1:6) .+ data.σs[1][1],
    color = my_black,
    linewidth = 4,
    linestyle = :dash,
)
Colorbar(fig[1, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
# ATTRACTIVE

data = load_object(
    "data/Thermal/Single_σ0[35]_σdot0[20]_MemInf_λ1_Φ-1_μ1_d60_ΩT1.0e-6_τ120.jld2",
)
δ = data.τs[2] - data.τs[1]
rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
mx = maximum(abs.(rr))
hm = heatmap!(
    ax2,
    data.τs[1:step_size:end],
    collect(1:size(data.ρs)[1]) .* data.α,
    rr[1:step_size:end, :],
    colormap = :RdBu,
    colorrange = (-mx, mx),
)
lines!(ax2, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)
# xlims!(ax2, (0, 80))
ylims!(ax2, (0, 2000))

lines!(
    ax2,
    0:0.1:6,
    π * 10 * 9 * (0:0.1:6) .+ data.σs[1][1],
    color = my_black,
    linewidth = 4,
    linestyle = :dash,
)
Colorbar(fig[2, 2], hm; label = L"\Delta\rho", width = 15, ticksize = 15, tickalign = 1)
fig
# save("General_Example.pdf", fig)

# δ = system.δ
# τ = 1.1 * (α / σ_dot)
# n_pts = τ / δ |> floor |> Int
# ρHs = zeros(nChain, n_pts)
# tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

# # res = motion_solver(system, Φ0, λ, α, [40], [50], μ, tTraj, Inf, τ)


# lines(res.τs, vec(res.σs))

# Φ0 = 1 / 2
# σ0 = 45
# nChain = 10

# Δs = zeros(vPts)


# @showprogress for ii = 1:vPts
#     σ_dot = σ_dots[ii]
#     τ = 1.1 * (α / σ_dot)
#     n_pts = τ / δ |> floor |> Int
#     ρHs = zeros(nChain, n_pts)
#     tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)

#     res = motion_solver(system, Φ0, 1 / 4, α, [σ0], [σ_dot], 1, tTraj, Inf, τ)
#     σs = res.σs |> vec
#     mid_pt_idx = findmin(abs.(σs .- (σ0 + α)))[2]
#     v_final = (σs[mid_pt_idx] - σs[mid_pt_idx-1]) / δ
#     Δs[ii] = μ / 2 / (2 * π)^2 * (σ_dot^2 - v_final^2)
# end

# λ = 1 / 4
# Φ = 1 / 2

# analytic = en_loss.(σ_dots)



# lines!(ax1,log.(σ_dots), log.(Δs))
# # lines!(ax1,log.(σ_dots), log.(Δs_2))
# lines!(ax1, log.(σ_dots), log.(analytic))
# Δs- Δs_2


# lines!(ax1, (σ_dots), (Δs))
# lines!(ax1,(σ_dots), (Δs_2))
# lines!(ax1, (σ_dots), (analytic))
fig


# ## ENERGY LOSS
# data = load_object(
#     "data/non_thermal/Single_σ0[60]_σdot0[25]_MemInf_λ1.0_Φ-1.0_μ1_d60_ΩTnothing_τ100.jld2",
# )
# λ = data.λ
# Ω = data.ωmax
# Φ = data.Φ
# λ = 1/1
# function en_loss(v)
#     z = (2 * π * λ / v)^2
#     return (
#         4 * π^3 * Φ^2 / v^2 *
#         z *
#         exp(-z * (Ω^2 + 1) / 2) *
#         (
#             besseli(0, z * (Ω^2 - 1) / 2) +
#             (Ω^2 - 1) / 2 * (besseli(0, z * (Ω^2 - 1) / 2) - besseli(1, z * (Ω^2 - 1) / 2))
#         )
#     )
# end
# en_loss(25)
# vs = range(4, 20, length = 1000)
# r = en_loss.(vs)
# lines(vs, r)

# σ = [x[1] for x in data.σs] |> vec
# σ_dot = (σ[2:end] - σ[1:end-1]) ./ (data.τs[2] - data.τs[1])
# idx = findall(x -> σ_dot[x] < σ_dot[x-1] && σ_dot[x] < σ_dot[x+1], 2:(length(σ_dot)-1))
# kin_en = σ_dot[idx] .^ 2 / 2 / (2 * π)^2
# Δ_kin_en = kin_en[2:end] - kin_en[1:end-1]
# lines(1:length(σ_dot), σ_dot)
# fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
# mean_speed = [mean(σ_dot[idx[i]:idx[i+1]]) for i in 1:length(idx)-1]
# scatter!(ax1,σ_dot[idx[2:end]], Δ_kin_en)
# scatter!(ax1, σ_dot[idx[2:end]], -en_loss.(σ_dot[idx[2:end]]))

# scatter!(ax1,mean_speed, Δ_kin_en)
# scatter!(ax1,mean_speed[1:70],  -en_loss.(mean_speed[1:70]))
# # σ_dot[250]
# # ylims!(ax1, -.025, 0)
# # σ_dot[idx[2:end]] - mean_speed
# # ( 25 .^ 2 / 2 / (2 * π)^2)-(24.908.^ 2 / 2 / (2 * π)^2)
# # en_loss.(mean_speed)
# xlims!(ax1, (10,25))
# # ylims!(ax)
# fig
# mean_speed[1:70]|>minimum
# σ_dot[1:100]
# scatter(1:400, σ_dot[1:400])

# lines(σ_dot[idx[1:end-1]], Δ_kin_en)


# lines(data.τs[idx], kin_en)

# # σ_dot_max = find


# lines(data.τs[2:end], σ_dot)




# # # POTENTIAL SHAPE

# # fig = Figure(resolution = (1200, 1600), font = "CMU Serif", fontsize = 36)
# # ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
# # ax2 = Axis(fig[2, 1], xlabel = L"\tau", ylabel = L"\sigma")
# # # ax3 = Axis(fig[3, 1], xlabel = L"\tau", ylabel = L"\dot{\sigma}")

# using SpecialFunctions
# # # REPULSIVE

# @time besseli(0, 1)
# data = load_object(
#     "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ1.0_Φ1.0_μ1_d60_ΩTnothing_τ100.jld2",
# )
# s = [x[1] for x in data.σs] |> vec
# lines(
#     (data.τs[2:end]),
#     log.(-((-s[2:end] + s[1:end-1]) .^ 2) .^ (2 / 10) .+ ((s[2] .- s[1])^2)^(2 / 10)),
#     color = my_black,
#     linewidth = 1,
# )
# # lines((data.τs[2:end]), (s[2:end] - s[1:end-1])./(data.τs[2] - data.τs[1]), color = my_black, linewidth = 1)
# # lines(data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)

# a
# # δ = data.τs[2] - data.τs[1]
# # rr = reduce(hcat, [data.ρs[ii, :] .- ii * data.α for ii = 1:size(data.ρs)[1]])
# # mx = maximum(abs.(rr))
# # hm = heatmap!(
# #     ax1,
# #     data.τs[1:step_size:end],
# #     collect(1:size(data.ρs)[1]) .* data.α,
# #     rr[1:step_size:end, :],
# #     colormap = :vik,
# #     colorrange = (-mx, mx),
# # )
# # lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)
# # # xlims!(ax1, (0, 80))
# # ylims!(ax1, (0, 1500))

# # # lines!(
# # #     ax1,
# # #     0:0.1:6,
# # #     π * 10 * 9 * (0:0.1:6) .+ data.σs[1][1],
# # #     color = my_black,
# # #     linewidth = 4,
# # #     linestyle = :dash,
# # # )
# # fig

# # a
# # # xlims!(ax2, (0, 80))
# # # # ylims!(ax3, (0, 1000))
# # # # vv[1:step_size:end, :]
# # # lines!(ax2,data.τs[2: end], (data.σs[2:end] .- data.σs[1:end-1]) ./ δ, linewidth = 4, color = my_vermillion)
# # # vv[1:step_size:end, :]

# # # data = load_object(
# # #     "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ1_Φ1_μ1_d60_ΩTnothing_τ80.jld2",
# # # )


# # # lines(
# # #     (abs.(data.τs[1:45000-1])),
# # #     (abs.(data.σs[2:45000] - data.σs[1:45000-1] |> vec)) ./ (data.τs[2] - data.τs[1]),
# # # )

# # # xs = range(0, 100, length = 100)
# # # ys = range(0, 100, length = 48000)
# # # rr = reduce(hcat, [data.ρs[ii, :] .- ii * 10 for ii = 1:100])


# # # xs = range(0, 10, length = 25)
# # # ys = range(0, 15, length = 35)
# # # zs = [cos(x) * sin(y) for x in xs, y in ys]

# # # heatmap(xs, ys, zs)

# # # a

# # # ## s dependence

# # # fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 18)

# # # ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
# # # ax2 = Axis(fig[1, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

# # # ax3 = Axis(fig[2, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
# # # ax4 = Axis(fig[2, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

# # # ax = [ax1, ax2, ax3, ax4]
# # # lims = [[(0, 15), (0, 30)], [(0, 15), (0, 30)], [(0, 15), (0, 30)], [(0, 15), (0, 30)]]
# # # files = [
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # # ]
# # # txt = ["(a)", "(b)", "(c)", "(d)"]
# # # for ii = 1:length(ax)
# # #     mkFigure(ax[ii], files[ii], lims[ii])
# # #     text!(ax[ii], txt[ii], position = (1/2, 27), textsize = 18, font = "CMU Serif")
# # # end

# # # fig
# # # save("s_Dependence.png", fig)


# # # fig = Figure(resolution = (600, 800), font = "CMU Serif", fontsize = 18)

# # # ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
# # # ax2 = Axis(fig[2, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

# # # ax = [ax1, ax2]
# # # lims = [[(0, 8), (0, 50)], [(0, 8), (0, 50)], [(0, 8), (0, 50)], [(0, 8), (0, 50)]]
# # # files = [
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[6]_MemInf_s0.125_F0.25_MInf_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[4]_MemInf_s0.125_F0.25_MInf_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # # ]
# # # txt = ["(a)", "(b)", "(c)", "(d)"]
# # # for ii = 1:length(ax)
# # #     mkFigure(ax[ii], files[ii], lims[ii])
# # #     text!(ax[ii], txt[ii], position = (1/2, 27), textsize = 18, font = "CMU Serif")
# # # end

# # # fig
# # # save("s_Dependence.png", fig)



# # # ## Speed dependence

# # # fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 18)

# # # ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
# # # ax2 = Axis(fig[1, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

# # # ax3 = Axis(fig[2, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
# # # ax4 = Axis(fig[2, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

# # # ax = [ax1, ax2, ax3, ax4]
# # # txt = ["(a)", "(b)", "(c)", "(d)"]

# # # files = [
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[1.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[2.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[4.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[5.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
# # # ]

# # # for ii = 1:length(ax)
# # #     mkFigure(ax[ii], files[ii], txt[ii])
# # # end

# # # fig
# # # ylims!(ax3, (0, 40))
# # # ylims!(ax4, (0, 40))
# # # save("s_F_Dependence.png", fig)




















# # # xlims!(ax2, (0, 20))

# # # fig
# # # [x[1] for x in data.Rs]







# # # fig = Figure(resolution = (2400, 3200), font = "CMU Serif", fontsize = 32)

# # # ax1 = Axis(fig[1, 1])
# # # ax2 = Axis(fig[1, 2])
# # # # ax3 = Axis(fig[1, 3])

# # # ax4 = Axis(fig[2, 1])
# # # ax5 = Axis(fig[2, 2])
# # # # ax6 = Axis(fig[2, 3])

# # # ax7 = Axis(fig[3, 1])
# # # ax8 = Axis(fig[3, 2])
# # # # ax9 = Axis(fig[3, 3])

# # # ax10 = Axis(fig[4, 1])
# # # ax11 = Axis(fig[4, 2])
# # # # ax12 = Axis(fig[4, 3])

# # # # ax13 = Axis(fig[5, 1])
# # # # ax14 = Axis(fig[5, 2])
# # # # ax15 = Axis(fig[5, 3])

# # # # ax16 = Axis(fig[6, 1])
# # # # ax17 = Axis(fig[6, 2])
# # # # ax18 = Axis(fig[6, 3])
# # # ax = [
# # #     ax1,
# # #     ax2,
# # #     # ax3,
# # #     ax4,
# # #     ax5,
# # #     # ax6,
# # #     ax7,
# # #     ax8,
# # #     # ax9,
# # #     ax10,
# # #     ax11,
# # #     # ax12,
# # #     # ax13,
# # #     # ax14,
# # #     # ax15,
# # #     # ax16,
# # #     # ax17,
# # #     # ax18,
# # # ]
# # # files = [
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.25_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.25_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.25_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.125_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.125_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.125_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F-0.125_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.125_M1_d60_ΩTnothing_τ20.jld2",
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.125_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F-0.25_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.25_M1_d60_ΩTnothing_τ20.jld2",
# # #     # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.25_M1_d60_ΩTnothing_τ20.jld2",
# # # ]

# # # for ii = 1:length(ax)
# # #     mkFigure(ax[ii], files[ii])

# # # end

# # # fig
# # # save("test.pdf", fig)



# # # ## Single Trajectory Analysis
# # # data = load_object(
# # #     "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.03125_M1_d60_ΩTnothing_τ20.jld2",
# # # )

# # # step = 50


# # # # l2 = lines!(ax, (1:length(Rs)) .* δ, [x[2] for x in res.Rs] .- 0, color = my_red)
# # # # l2 = lines!(ax, (1:length(Rs)) .* δ, [x[3] for x in res.Rs] .- 0, color = my_green)
# # # fig

# # # # data.rs[1]



# # # m = data.m                      # Mass of the chain atoms
# # # k = data.k                      # Spring force constant
# # # K = data.K                      # Confining potential force constant
# # # Ωmin = sqrt(K / m)              # Smallest chain frequency
# # # Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency

# # # M = data.M                      # Mass of the mobile particles
# # # K_M = data.K_M                  # Trap force constant
# # # ΩM = √(K_M / M)                 # Trap frequency
# # # t_M = 2 * π / ΩM                # Period of the trapped mass
# # # Rs = data.Rs |> vec
# # # rs = data.rs
# # # ts = data.ts ./ t_M

# # # let
# # #     fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

# # #     ax1 = Axis(fig[1, 1], xlabel = L"t/t_M", ylabel = L"R")
# # #     ax3 = Axis(fig[1, 2], xlabel = L"t/t_M", ylabel = L"R")
# # #     ax5 = Axis(fig[1, 3], xlabel = L"t/t_M", ylabel = L"R")

# # #     ax2 = Axis(fig[2, 1], xlabel = L"t/t_M", ylabel = L"r")
# # #     ax4 = Axis(fig[2, 2], xlabel = L"t/t_M", ylabel = L"r")
# # #     ax6 = Axis(fig[2, 3], xlabel = L"t/t_M", ylabel = L"r")

# # #     ax_x_lims = [(0, 120), (0, 120), (88, 92), (88, 92), (116, 117), (116, 117)]
# # #     ax_y_lims =
# # #         [(-11, 11), (-0.15, 0.15), (-5.5, 5.5), (-0.15, 0.15), (-0.04, 0.04), (-0.04, 0.04)]
# # #     ax = [ax1, ax2, ax3, ax4, ax5, ax6]
# # #     txt = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

# # #     for n = 1:length(ax)
# # #         xlims!(ax[n], ax_x_lims[n])
# # #         ylims!(ax[n], ax_y_lims[n])
# # #         text!(
# # #             ax[n],
# # #             txt[n],
# # #             position = (
# # #                 ax_x_lims[n][1] + 5 / 6 * (ax_x_lims[n][2] - ax_x_lims[n][1]),
# # #                 ax_y_lims[n][2] / 2,
# # #             ),
# # #             textsize = 16,
# # #             font = "CMU Serif",
# # #         )
# # #     end
# # #     idx = findall(x -> (x <= ax_x_lims[1][2]) && (x >= ax_x_lims[1][1]), ts)
# # #     lines!(ax[1], ts[idx], Rs[idx], color = my_red, linewidth = 1)
# # #     lines!(ax[2], ts[idx], rs[idx], color = my_red, linewidth = 1)

# # #     idx = findall(x -> (x <= ax_x_lims[3][2]) && (x >= ax_x_lims[3][1]), ts)
# # #     lines!(ax[3], ts[idx], Rs[idx], color = my_red, linewidth = 1)
# # #     lines!(ax[4], ts[idx], rs[idx], color = my_red, linewidth = 1)


# # #     idx = findall(x -> (x <= ax_x_lims[5][2]) && (x >= ax_x_lims[5][1]), ts)
# # #     lines!(ax[5], ts[idx], Rs[idx], color = my_red, linewidth = 1)
# # #     lines!(ax[6], ts[idx], rs[idx], color = my_red, linewidth = 1)

# # #     colgap!(fig.layout, 5)
# # #     rowgap!(fig.layout, 5)
# # #     save("General_Example.pdf", fig)
# # # end
