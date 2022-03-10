include("../src/main.jl")

colors = [my_red, my_green, my_blue, my_violet]
### keeping F fixed at 0.25
s = [1 / 16, 1 / 8, 1 / 4, 1 / 2]
data = [
    load_object(
        "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.0625_F-0.25_M0.5_d60_ΩTnothing_τ100.jld2",
    ),
    load_object(
        "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.125_F-0.25_M0.5_d60_ΩTnothing_τ100.jld2",
    ),
    load_object(
        "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.25_F-0.25_M0.5_d60_ΩTnothing_τ100.jld2",
    ),
    load_object(
        "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.5_F-0.25_M0.5_d60_ΩTnothing_τ100.jld2",
    ),
]

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 14)
ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
ax3 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
ax4 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
axes = [ax1, ax2, ax3, ax4]
colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)


#check for width scaling (potential width)
labs = [L"s = 1/16", L"s = 1/8", L"s = 1/4", L"s = 1/2"]

Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
tss = [data[i].ts ./ tmin for i = 1:length(data)]


for ii = 1:length(axes)
    lines!(axes[ii], tss[ii], [x[1] for x in data[ii].Rs])
    text!(axes[ii], "s = $(s[ii])", position = (40, 2), textsize = 30, font = "CMU Serif")
end
supertitle = Label(fig[0, :], "F = -0.25", textsize = 30)
fig

F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]
s = [1 / 16, 1 / 8, 1 / 4, 1 / 2]
for f in F
    data = [
        load_object(
            "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s$(z)_F$(f)_M0.5_d60_ΩTnothing_τ100.jld2",
        ) for z in s
    ]
    println(length(data))
    Ωmin = √(data[1].K / data[1].m)       # Minimum phonon frequency
    tmin = 2 * π / Ωmin
    tss = [data[1].ts ./ tmin]           # time steps should be same for all
    fig = Figure(resolution = (2000, 1200), font = "CMU Serif", fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
    ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
    ax3 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
    ax4 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
    axes = [ax1, ax2, ax3, ax4]
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)

    for ii = 1:length(data)
        lines!(axes[ii], tss[1], [x[1] for x in data[ii].Rs])
        text!(
            axes[ii],
            "s = $(s[ii])",
            position = (80, 5),
            textsize = 30,
            font = "CMU Serif",
        )
    end
    supertitle = Label(fig[0, :], "F = ($f)", textsize = 30)
    display(fig)
end


#check for depth dependence
F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]
s = [1 / 16, 1 / 8, 1 / 4, 1 / 2]
for z in s
    data = [
        load_object(
            "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s$(z)_F$(f)_M0.5_d60_ΩTnothing_τ100.jld2",
        ) for f in F
    ]
    Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
    tmin = 2 * π / Ωmin
    tss = [data[i].ts ./ tmin for i = 1:length(data)]
    fig = Figure(resolution = (2000, 1200), font = "CMU Serif", fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
    ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
    ax3 = Axis(fig[1, 3], xlabel = L"ts / t_min", ylabel = L"R")
    ax4 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
    ax5 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
    ax6 = Axis(fig[2, 3], xlabel = L"ts / t_min", ylabel = L"R")
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)

    for ii = 1:length(data)
        lines!(axes[ii], tss[ii], [x[1] for x in data[ii].Rs])
        text!(
            axes[ii],
            "F = $(F[ii])",
            position = (80, 5),
            textsize = 30,
            font = "CMU Serif",
        )
    end
    supertitle = Label(fig[0, :], "s = ($z)", textsize = 30)
    display(fig)
end


###s = 0.0625
# F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]

# data = [
#     load_object(
#         "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.0625_F$(f)_M0.5_d60_ΩTnothing_τ100.jld2",
#     ) for f in F
# ]

# Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
# tmin = 2 * π / Ωmin                 # Period of the slowest mode
# tss = [data[i].ts ./ tmin for i = 1:length(data)]
# fig = Figure(resolution = (2000, 1200), font = "CMU Serif", fontsize = 20)
# ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax3 = Axis(fig[1, 3], xlabel = L"ts / t_min", ylabel = L"R")
# ax4 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax5 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax6 = Axis(fig[2, 3], xlabel = L"ts / t_min", ylabel = L"R")
# axes = [ax1, ax2, ax3, ax4, ax5, ax6]
# colgap!(fig.layout, 5)
# rowgap!(fig.layout, 5)

# for ii = 1:length(axes)
#     lines!(axes[ii], tss[ii], [x[1] for x in data[ii].Rs])
#     text!(axes[ii], "F = $(F[ii])", position = (80, 5), textsize = 30, font = "CMU Serif")
# end

# fig

# ### s = 0.125

# F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]

# data = [
#     load_object(
#         "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.125_F$(f)_M0.5_d60_ΩTnothing_τ100.jld2",
#     ) for f in F
# ]

# Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
# tmin = 2 * π / Ωmin                 # Period of the slowest mode
# tss = [data[i].ts ./ tmin for i = 1:length(data)]
# fig = Figure(resolution = (2000, 1200), font = "CMU Serif", fontsize = 20)
# ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax3 = Axis(fig[1, 3], xlabel = L"ts / t_min", ylabel = L"R")
# ax4 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax5 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax6 = Axis(fig[2, 3], xlabel = L"ts / t_min", ylabel = L"R")
# axes = [ax1, ax2, ax3, ax4, ax5, ax6]
# colgap!(fig.layout, 5)
# rowgap!(fig.layout, 5)

# for ii = 1:length(axes)
#     lines!(axes[ii], tss[ii], [x[1] for x in data[ii].Rs])
#     text!(axes[ii], "F = $(F[ii])", position = (80, 5), textsize = 30, font = "CMU Serif")
# end
# supertitle = Label(fig[0, :], "s = 0.125", textsize = 30)
# fig

# ### s = 0.25

# F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]

# data = [
#     load_object(
#         "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.25_F$(f)_M0.5_d60_ΩTnothing_τ100.jld2",
#     ) for f in F
# ]

# Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
# tmin = 2 * π / Ωmin                 # Period of the slowest mode
# tss = [data[i].ts ./ tmin for i = 1:length(data)]
# fig = Figure(resolution = (2000, 1200), font = "CMU Serif", fontsize = 20)
# ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax3 = Axis(fig[1, 3], xlabel = L"ts / t_min", ylabel = L"R")
# ax4 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax5 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax6 = Axis(fig[2, 3], xlabel = L"ts / t_min", ylabel = L"R")
# axes = [ax1, ax2, ax3, ax4, ax5, ax6]
# colgap!(fig.layout, 5)
# rowgap!(fig.layout, 5)

# for ii = 1:length(axes)
#     lines!(axes[ii], tss[ii], [x[1] for x in data[ii].Rs])
#     text!(axes[ii], "F = $(F[ii])", position = (80, 5), textsize = 30, font = "CMU Serif")
# end
# supertitle = Label(fig[0, :], "s = 0.25", textsize = 30)
# fig

# ### s = 0.5

# F = [-1, -1 / 2, -1 / 4, 1 / 4, 1 / 2, 1]

# data = [
#     load_object(
#         "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.5_F$(f)_M0.5_d60_ΩTnothing_τ100.jld2",
#     ) for f in F
# ]

# Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
# tmin = 2 * π / Ωmin                 # Period of the slowest mode
# tss = [data[i].ts ./ tmin for i = 1:length(data)]
# fig = Figure(resolution = (2000, 1200), font = "CMU Serif", fontsize = 20)
# ax1 = Axis(fig[1, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax2 = Axis(fig[1, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax3 = Axis(fig[1, 3], xlabel = L"ts / t_min", ylabel = L"R")
# ax4 = Axis(fig[2, 1], xlabel = L"ts / t_min", ylabel = L"R")
# ax5 = Axis(fig[2, 2], xlabel = L"ts / t_min", ylabel = L"R")
# ax6 = Axis(fig[2, 3], xlabel = L"ts / t_min", ylabel = L"R")
# axes = [ax1, ax2, ax3, ax4, ax5, ax6]
# colgap!(fig.layout, 5)
# rowgap!(fig.layout, 5)

# for ii = 1:length(axes)
#     lines!(axes[ii], tss[ii], [x[1] for x in data[ii].Rs])
#     text!(axes[ii], "F = $(F[ii])", position = (80, 5), textsize = 30, font = "CMU Serif")
# end
# supertitle = Label(fig[0, :], "s = 0.5", textsize = 30)
# hidexdecorations!(ax1, ticks = false)
# #hideydecorations!(ax1, ticks = false)
# fig


# plot for mass dependence
F = [1.0]
s = [1 / 16]
M = [0.25, 0.5, 1, 1.25, 1.5, 2.0]
data = [
    load_object(
        "../data/Non_Thermal/Single_x0[5]_v0[0.5]_MemInf_s0.0625_F1.0_M$(w)_d60_ΩTnothing_τ100.jld2",
    ) for w in M
]
Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
tss = data[1].ts ./ tmin
fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 20)
ax1 = Axis(fig[1, 1], xlabel = L"t / t_{min}", ylabel = L"R")
ax2 = Axis(fig[1, 2], xlabel = L"t / t_{min}", ylabel = L"R")
ax3 = Axis(fig[1, 3], xlabel = L"t / t_{min}", ylabel = L"R")
ax4 = Axis(fig[2, 1], xlabel = L"t / t_{min}", ylabel = L"R")
ax5 = Axis(fig[2, 2], xlabel = L"t / t_{min}", ylabel = L"R")
ax6 = Axis(fig[2, 3], xlabel = L"t / t_{min}", ylabel = L"R")
axes = [ax1, ax2, ax3, ax4, ax5, ax6]
for ii=1:length(data)
    lines!(axes[ii], tss, [x[1] for x in data[ii].Rs])
    text!(axes[ii], "M = $(M[ii])", position = (5,5), textsize = 30, font = "CMU Serif")
end
supertitle = Label(fig[0, :], "s = $(s[1]), F = $(F[1])", textsize = 30)
fig


### plot for velocity dependence
F = [1.0]
s = [1/16]
M = [1.0]
v0 = v0 = [[1 / 4], [1 /2] , [1.0], [3 / 2], [2.0]]
data = [
    load_object(
        "../data/Non_Thermal/Single_x0[5]_v0$(v)_MemInf_s0.0625_F1.0_M1.0_d60_ΩTnothing_τ100.jld2",
    ) for v in v0
]
Ωmin = √(data[4].K / data[4].m)       # Minimum phonon frequency
tmin = 2 * π / Ωmin                 # Period of the slowest mode
tss = data[1].ts ./ tmin
fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 20)
ax1 = Axis(fig[1, 1], xlabel = L"t / t_{min}", ylabel = L"R")
ax2 = Axis(fig[1, 2], xlabel = L"t / t_{min}", ylabel = L"R")
ax3 = Axis(fig[1, 3], xlabel = L"t / t_{min}", ylabel = L"R")
ax4 = Axis(fig[2, 1], xlabel = L"t / t_{min}", ylabel = L"R")
ax5 = Axis(fig[2, 2], xlabel = L"t / t_{min}", ylabel = L"R")
axes = [ax1, ax2, ax3, ax4, ax5, ax6]
for ii=1:length(data)
    lines!(axes[ii], tss, [x[1] for x in data[ii].Rs])
    text!(axes[ii], "v = $(v0[ii])", position = (5,5), textsize = 30, font = "CMU Serif")
end
supertitle = Label(fig[0, :], "s = $(s[1]), F = $(F[1])", textsize = 30)
fig


### rerun for s = 1/16, F = 1/4, m = 0.5, v = 0.5
data = load_object("../data/Non_Thermal/rerun_Single_x0[5]_v[0.5]_MemInf_s0.0625_F0.25_M0.5_d60_ΩTnothing_τ100.jld2")
Ωmin = √(data.K / data.m)
tmin = 2 * π / Ωmin 
tss = data.ts ./ tmin
fig = Figure(resolution = (800,800), font = "CMU Serif", fontsize = 20)
ax1 = Axis(fig[1,1], xlabel = L"t / t_{min}", ylabel = "R")
lines!(ax1, tss, [x[1] for x in data.Rs])
supertitle = Label(fig[0,:], "s = 0.0625, F = 0.5, M = 0.5, v = 0.5")
fig
