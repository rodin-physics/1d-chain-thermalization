include("../src/main.jl")
# @everywhere include("../src/main.jl")

system = load_object("precomputed/systems/System_ωmax10_d60_l200.jld2")

d = 60
τ = 80                            # Simulation time
δ = system.δ                        # Time step
a = 20                               # Distance between chain atoms
n_pts = τ / δ |> floor |> Int
nChain = 100
ρHs = zeros(nChain, n_pts)
tTraj = ThermalTrajectory(system.ωmax, system.δ, ρHs, nothing)
λ = 1;
Φ = 1;
@time res = motion_solver(system, 1, 1, a, [25], [20], 1, tTraj, 2, τ)
@time res_N = motion_solver(system, 1.01, 3, a, [25], [20], 1, tTraj, 2, τ)
@time res_NN = motion_solver(system, 1.965, 6, a, [25], [20], 1, tTraj, 2, τ)
@time res_NNN = motion_solver(system, 12, 9, a, [25], [20], 1, tTraj, 2, τ)
@time res_NNNN = motion_solver(system, 203, 12, a, [25], [20], 1, tTraj, 2, τ)



@time res_S = motion_solver(system, 1, 1, a, [res.σs[23450]], [(res.σs[23450]-res.σs[23449]) / (res.τs[2]-res.τs[1])], 1, tTraj, Inf, τ-res.τs[23450])
Plots.plot((res.τs), (res.σs |> vec))
Plots.plot!((res_N.τs), (res_N.σs |> vec))
Plots.plot!((res_NN.τs), (res_NN.σs |> vec))
Plots.plot!((res_NNN.τs), (res_NNN.σs |> vec))
Plots.plot!((res_NNNN.τs), (res_NNNN.σs |> vec))
Plots.plot!((res_S.τs).+res.τs[23450], (res_S.σs |> vec))

(res.σs[23450]-res.σs[23449]) / (res.τs[2]-res.τs[1])
res.τs[23450]
Plots.plot((res_S.τs), (res_S.σs |> vec))


Plots.plot!((res.τs), (res.σs |> vec))
Plots.plot((res.τs), (res.σs |> vec) -  (res_S.σs |> vec))

# @time res = motion_solver(system, Φ/5, λ/32, a, [25], [4], 1, tTraj, Inf, τ)
# @time res_test = motion_solver_test(system, Φ, λ, a, [25], [10], 1, tTraj, Inf, τ)
# println(res == res_test)
# motion_solver(system, Φ, λ, a, [25], [10], 1, tTraj, Inf, τ) ==
# motion_solver(system, Φ, λ, a, [25], [10], 1, tTraj, Inf, τ)
# res.ρs ≈ res_test.ρs
# res
# @time res = motion_solver(system, 10, 1, a, [25], [5], 1, tTraj, 50, τ)
1/2*(4 / 2 / pi)^2
# using Plots
# res.σs[1]
# res.σs[3]
# # using Plots
# # motion_solver_test()
# maximum(res.σs |> vec)
# res.σs[2]-res.σs[1]
40/600
using Plots

Plots.plot((res.τs), (res.ρs[10,:] |> vec))
length(res.τs)
a
# Plots.plot((abs.(res.τs[1:end-1])), (abs.(res.σs[2:end] - res.σs[1:end-1] |> vec)))




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