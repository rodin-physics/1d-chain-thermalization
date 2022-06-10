include("../src/main.jl")
# @everywhere include("../src/main.jl")

system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")

d = 60
τ = 80                            # Simulation time
δ = system.δ                        # Time step
a = 20                               # Distance between chain atoms
n_pts = τ / δ |> floor |> Int
nChain = 50
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
λ = 1;
Φ = 1;
@time res1 = motion_solver(system, 1, 1, a, [25], [5], 1, tTraj, Inf, τ)
# @time res1 = motion_solver(system, 1, 1, a, [25], [5], 1, tTraj, Inf, τ)
# @time res = motion_solver_TEST(system, 1, 1, a, [25], [5], 1, tTraj, Inf, τ)
@time res = motion_solver_TEST(system, 1, 1, a, [25], [5], 1, tTraj, Inf, τ)
# println(res1.σs ≈res.σs)
# res.σs - res1.σs
# a
# @time res_N = motion_solver(system, 1.01, 3, a, [25], [20], 1, tTraj, 2, τ)
# @time res_NN = motion_solver(system, 1.965, 6, a, [25], [20], 1, tTraj, 2, τ)
# @time res_NNN = motion_solver(system, 12, 9, a, [25], [20], 1, tTraj, 2, τ)
# @time res_NNNN = motion_solver(system, 203, 12, a, [25], [20], 1, tTraj, 2, τ)



# system = load_object("precomputed/systems/System_ωmax10_d3000_l10.jld2")

# d = 3000
# τ = 0.45                            # Simulation time
# δ = system.δ                        # Time step
# a = 10                               # Distance between chain atoms
# n_pts = τ / δ |> floor |> Int
# nChain = 10
# ρHs = zeros(nChain, n_pts)
# tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
# λ = 1;
# Φ = -1;
# μ = 1;
# @time res = motion_solver(system, Φ, λ, a, [30], [25], μ, tTraj, Inf, τ)
# s = (res.σs |> vec)
# v = (s[2:end] - s[1:end-1]) / (res.τs[2] - res.τs[1])
# scatter((res.τs)[2:end], v)






# kin_en = μ * v .^ 2 / 2 / (2 * π)^2
# kin_en[end] - kin_en[1]

# Δ_kin_en = kin_en[2:end] - kin_en[1:end-1]


# @time res_S = motion_solver(
#     system,
#     1,
#     1,
#     a,
#     [res.σs[23450]],
#     [(res.σs[23450] - res.σs[23449]) / (res.τs[2] - res.τs[1])],
#     1,
#     tTraj,
#     Inf,
#     τ - res.τs[23450],
# )
# Plots.plot((res.τs), (res.σs |> vec))
# Plots.plot!((res_N.τs), (res_N.σs |> vec))
# Plots.plot!((res_NN.τs), (res_NN.σs |> vec))
# Plots.plot!((res_NNN.τs), (res_NNN.σs |> vec))
# Plots.plot!((res_NNNN.τs), (res_NNNN.σs |> vec))
# Plots.plot!((res_S.τs) .+ res.τs[23450], (res_S.σs |> vec))

# (res.σs[23450] - res.σs[23449]) / (res.τs[2] - res.τs[1])
# res.τs[23450]
# Plots.plot((res_S.τs), (res_S.σs |> vec))


# Plots.plot!((res.τs), (res.σs |> vec))
# Plots.plot((res.τs), (res.σs |> vec) - (res_S.σs |> vec))
# @time Γ(4, 0:10, 10);
# # @time res = motion_solver(system, Φ/5, λ/32, a, [25], [4], 1, tTraj, Inf, τ)
# # @time res_test = motion_solver_test(system, Φ, λ, a, [25], [10], 1, tTraj, Inf, τ)
# # println(res == res_test)
# # motion_solver(system, Φ, λ, a, [25], [10], 1, tTraj, Inf, τ) ==
# # motion_solver(system, Φ, λ, a, [25], [10], 1, tTraj, Inf, τ)
# # res.ρs ≈ res_test.ρs
# # res
# # @time res = motion_solver(system, 10, 1, a, [25], [5], 1, tTraj, 50, τ)
# 1 / 2 * (10 / 2 / pi)^2
# 1 / 2 * (20 / 2 / pi)^2
# 1 / 2 * (7 / 2 / pi)^2
# # using Plots
# # res.σs[1]
# # res.σs[3]
# # # using Plots
# # # motion_solver_test()
# # maximum(res.σs |> vec)
# # res.σs[2]-res.σs[1]
# 40 / 600
# using Plots

# Plots.plot((res.τs), (res.ρs[10, :] |> vec))
# length(res.τs)
# a
# # Plots.plot((abs.(res.τs[1:end-1])), (abs.(res.σs[2:end] - res.σs[1:end-1] |> vec)))




# scatter((res.τs[1:end-1]),(0.5*( (res.σs[2:end] - res.σs[1:end-1]) / 2 / pi / (res.τs[2] - res.τs[1])).^2 |> vec), markersize = 0.25)

# length(res.τs)
# lines(log.(abs.(res.τs[1:25000-1])), log.((res.σs[2] - res.σs[1]).^5 .-(res.σs[2:25000] - res.σs[1:25000-1] |> vec).^5))
# lines((abs.(res.τs[1:25000-1])), ((res.σs[2] - res.σs[1]).^(5/2) .-(res.σs[2:25000] - res.σs[1:25000-1] |> vec).^(5/2)))
# lines((abs.(data.τs[1:48000-1])), ((data.σs[2] - data.σs[1]).^(5/2) .-((data.σs[2:48000] - data.σs[1:48000-1] |> vec).^2).^(5/4)))
# ((res.σs[2] - res.σs[1]) / (res.τs[2]-res.τs[1]) )^2/2

# lines((abs.(res.τs[1:45000-1])), (abs.(res.σs[2:45000] - res.σs[1:45000-1] |> vec)))
# Plots.plot(res.τs, res.ρs[15,:] |> vec)
# lines(res.τs, res.ρs[19,:] |> vec)
# res.σs
# @time res = motion_solver(system, 1, 1, 1, [5.5], [0.25], 1, tTraj, Inf, τ)
# Base.summarysize(res)
# res.σs

# Plots.plot(1:50,res[:,2000])
# using Plots
# aa = ones(500000, 500);
# bb = ones(500000, 500);
# cc = ones(500000, 500);

# @time aa .= aa .+ aa;

# @time for ii = 1:500
#     Threads.@spawn(view(bb, 1:100000, ii) .= view(bb, 1:100000, ii) .+3 .* view(bb, 1:100000, ii))
# end

# @time Threads.@threads for ii = 1:500
#     view(cc, 1:100000, ii) .= view(cc, 1:100000, ii) .+  view(cc, 1:100000, ii)
# end


# # @time aa .= aa .* 3;

# @time for ii in 1 : 500
#     Threads.@spawn(view(bb,:, ii) .= view(bb,:, ii) .* 3)
# end


# @time Threads.@threads for ii in 1 : 500
#     view(cc,:, ii) .=  view(cc,:, ii) .* 3
# end



# proc_num = 9
# addprocs(proc_num - nprocs())



# R = Γ.(range(0,1/2,length = 1000), [0], 10)
# lines(range(0,1/2,length = 1000), R)

# fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 36)
# ax1 = Axis(fig[1, 1], xlabel = L"\tau", ylabel = L"\sigma")
# data = load_object(
#     "data/non_thermal/Single_σ0[55]_σdot0[20]_MemInf_λ1.0_Φ1.0_μ1_d60_ΩTnothing_τ100.jld2",
# )
# δ = data.τs[2] - data.τs[1]
# # lines!(ax1, data.τs .+0, ([x[1] for x in data.σs] |> vec) .+0, color = my_red, linewidth = 5)
# # lines!(ax1,data.τs[2: end] .+22.0, (data.σs[2:end] .- data.σs[1:end-1]) ./ δ, linewidth = 1, color = my_vermillion)
# v = (data.σs[2:end] .- data.σs[1:end-1]) ./ δ
# idx = findall(i -> v[i] > v[i-1] && v[i] > v[i+1], 2:length(v)-2)
# ch = v[idx[2:end]] - v[idx[1:end-1]]
# lines!(ax1, data.τs[2:end], (v[1] .- v) .^ (1 / 2))
# scatter(1:(length(idx)-0), (v[1] .- v[idx]) .^ (1 / 1))
# scatter(((v[idx[2:end]])[1:end-23]), (-ch[1:end-23]))
# # length(v)
# # v[10000]
# aa = (v[2:end] .- v[1:end - 1]) ./ δ;
# scatter!(ax1,data.τs[3: end] ,(aa))
# xlims!(ax1, (0,5))
# # lines!(ax1,data.τs[2: end] .+10.0, (data.σs[2:end] .- data.σs[1:end-1]) ./ δ, linewidth = 1, color = my_vermillion)

# data = load_object(
#     "data/non_thermal/Single_σ0[60]_σdot0[25]_MemInf_λ0.5_Φ-1.0_μ1_d60_ΩTnothing_τ100.jld2",
# )
# δ = data.τs[2] - data.τs[1]
# # v = (data.σs[2:end] .- data.σs[1:end-1]) ./ δ
# # aa = (v[2:end] .- v[1:end - 1]) ./ δ;
# # scatter!(ax1,data.τs[3: end].-10 ,aa, color = my_green)
# # lines!(ax1, data.τs, [x[1] for x in data.σs] |> vec, color = my_black, linewidth = 5)

# # lines!(ax1, data.τs[2: end], (data.σs[2:end] .- data.σs[1:end-1]) ./ δ, linewidth = 1, color = my_blue)
# fig

# gg = system.Γ[1, :]

# lines(data.τs[1:400], gg[1:400])
# lines(data.τs[1:4000], g)


# function GG(τ, ωmax)
#     int_fun(x) = cos(2 * π * τ * ω(ωmax, x)) / ω(ωmax, x)^2
#     res = quadgk(int_fun, 0, π / 2)
#     return (res[1] * 2 / π)
# end

# g = GG.(data.τs[1:4000], 5)
ts = (-1/10):(1/600):(1/10)
p = [(t * 25) .* exp(- (t * 25)^2/ (2 * 0.5^2)) for t in ts]
scatter(ts .* 25, p)