include("../src/main.jl")

# τmax = 100                         # Simulation time
# ωmax = 10
# d = 60
# δ = (1 / ωmax) / d              # Time step
# n_pts = floor(τmax / δ) |> Int  # Number of time steps given t_max and δ
# n_masses = 10000                  # Number of chain masses for simulating ρ0

# qa_s = 2 * pi .* (1:n_masses) / n_masses
# ωs = ω.(ωmax, qa_s ./ 2)
# ωT = 1e-5

# ζs = ζq.(ωs, ωT)
# ϕs = 2 * π * rand(n_masses)
# @benchmark ρH(100, δ, ζs, ϕs, ωs, 1:200)|>real

# @time ρH(100, δ, ζs, ϕs, ωs, 1:200)|>real
# @time r1 = [ζH(n, δ, ζs, ϕs, ωs, 1:100) for n in 1:1000];
# @time ζH(100, δ, ζs, ϕs, ωs, 1:200)|>real
# .06*120000
# # function ρH(n, δ, ζs, ϕs, ωs, gs)
# #     n_ζ = length(ζs)
# #     res =
# #         [
# #             exp.(1im * 2 * pi / n_ζ * x * gs .-1im * 2 * π * δ * n * ωs[x] .- 1im * ϕs[x]) * ζs[x] / √(n_ζ) for
# #             x = 1:n_ζ
# #         ] |> sum
# #     return res
# # end

# # function ρH(n, δ, ζs, ϕs, ωs, gs)
# #     n_ζ = length(ζs)
# #     res =
# #         [
# #             exp.(1im * 2 * π / n_ζ * x * gs) * ζs[x] *
# #             exp(-1im * (2 * π * δ * n * ωs[x] + ϕs[x])) / √(n_ζ) for
# #             x = 1:n_ζ
# #         ] |> sum
# #     return res
# # end

# function ρH(n, δ, ζs, ϕs, ωs, gs)
#     n_ζ = length(ζs)
#     f = transpose(exp.(-1im * (2 * π * δ * n * ωs + ϕs)) / √(n_ζ) .* ζs)
#     r = [f * exp.(1im * 2 * π / n_ζ .* (1:n_ζ) * g) for g in gs]
#     return r
# end


# function ρH_2(n, δ, ζs, ϕs, ωs, gs)
#     n_ζ = length(ζs)
#     f = exp.(-1im * (2 * π * δ * n * ωs + ϕs)) / √(n_ζ) .* ζs
#     p = [dot(f, exp.(1im * 2 * π / n_ζ .* (1:n_ζ) * g)) for g in gs]

#     return p
#     # res =
#     #     [
#     #         exp.(1im * 2 * π / n_ζ * x * gs) * ζs[x] *
#     #         exp(-1im * (2 * π * δ * n * ωs[x] + ϕs[x])) / √(n_ζ) for
#     #         x = 1:n_ζ
#     #     ] |> sum
#     # return res
# end
# @benchmark s = ρH_2(1, δ, ζs, ϕs, ωs, 1:200)
# @time d = ρH(1, δ, ζs, ϕs, ωs, 1:200);
# s ≈ d
# s == d
# s
# d

# using BenchmarkTools