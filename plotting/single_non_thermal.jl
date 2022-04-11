include("../src/main.jl")

step_size = 10

function mkFigure(ax, filename, lms)
    data = load_object(filename)
    t_slow = 2 * π / √(data.K / data.m)

    for ii = 1:length(data.rs[1])
        lines!(
            ax,
            data.ts[1:step_size:end]./ t_slow,
            # [x[ii] for x in data.rs[1:step_size:end]],
            ii * ones(length(data.ts[1:step_size:end])),
            colormap = :balance,
            color = [x[ii] for x in data.rs[1:step_size:end]] .- ii,
            colorrange = (-8e-2, 8e-2),
            linewidth = 1,
        )
    end
    lines!(
        ax,
        data.ts ./ t_slow,
        [x[1] for x in data.Rs],
        color = colorant"rgba(0, 0, 0, 0.75)",
        linewidth = 2,
    )
    xlims!(ax, lms[1])
    ylims!(ax, lms[2])
end



## General Example

fig = Figure(resolution = (600, 400), font = "CMU Serif", fontsize = 18)

ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
mkFigure(
    ax1,
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.5_F-0.5_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.0]_v0[0.5]_MemInf_s0.5_F-1_M1.0_d60_ΩTnothing_τ50.jld2",
    
    [(0, 20), (0, 60)],
)

lines!(
    ax1,
    0:0.1:2,
    2 * pi * 4 .* (0:0.1:2) .+ 5.5,
    color = colorant"rgba(0, 0, 0, 0.75)",
    linewidth = 2,
    linestyle = :dash,
)
fig

save("General_Example.pdf", fig)

## s dependence

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 18)

ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
ax2 = Axis(fig[1, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

ax3 = Axis(fig[2, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
ax4 = Axis(fig[2, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

ax = [ax1, ax2, ax3, ax4]
lims = [[(0, 15), (0, 30)], [(0, 15), (0, 30)], [(0, 15), (0, 30)], [(0, 15), (0, 30)]]
files = [
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
]
txt = ["(a)", "(b)", "(c)", "(d)"]
for ii = 1:length(ax)
    mkFigure(ax[ii], files[ii], lims[ii])
    text!(ax[ii], txt[ii], position = (1/2, 27), textsize = 18, font = "CMU Serif")
end

fig
save("s_Dependence.png", fig)


fig = Figure(resolution = (600, 800), font = "CMU Serif", fontsize = 18)

ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
ax2 = Axis(fig[2, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

ax = [ax1, ax2]
lims = [[(0, 8), (0, 50)], [(0, 8), (0, 50)], [(0, 8), (0, 50)], [(0, 8), (0, 50)]]
files = [
    "data/Non_Thermal/Single_x0[5.5]_v0[6]_MemInf_s0.125_F0.25_MInf_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[4]_MemInf_s0.125_F0.25_MInf_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
]
txt = ["(a)", "(b)", "(c)", "(d)"]
for ii = 1:length(ax)
    mkFigure(ax[ii], files[ii], lims[ii])
    text!(ax[ii], txt[ii], position = (1/2, 27), textsize = 18, font = "CMU Serif")
end

fig
save("s_Dependence.png", fig)



## Speed dependence

fig = Figure(resolution = (1200, 800), font = "CMU Serif", fontsize = 18)

ax1 = Axis(fig[1, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
ax2 = Axis(fig[1, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

ax3 = Axis(fig[2, 1], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")
ax4 = Axis(fig[2, 2], xlabel = L"t/t_\mathrm{Slow}", ylabel = L"x")

ax = [ax1, ax2, ax3, ax4]
txt = ["(a)", "(b)", "(c)", "(d)"]

files = [
    "data/Non_Thermal/Single_x0[5.5]_v0[1.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[2.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[4.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[5.0]_MemInf_s0.125_F1.0_M1_d60_ΩTnothing_τ20.jld2",
]

for ii = 1:length(ax)
    mkFigure(ax[ii], files[ii], txt[ii])
end

fig
ylims!(ax3, (0, 40))
ylims!(ax4, (0, 40))
save("s_F_Dependence.png", fig)




















xlims!(ax2, (0, 20))

fig
[x[1] for x in data.Rs]







fig = Figure(resolution = (2400, 3200), font = "CMU Serif", fontsize = 32)

ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
# ax3 = Axis(fig[1, 3])

ax4 = Axis(fig[2, 1])
ax5 = Axis(fig[2, 2])
# ax6 = Axis(fig[2, 3])

ax7 = Axis(fig[3, 1])
ax8 = Axis(fig[3, 2])
# ax9 = Axis(fig[3, 3])

ax10 = Axis(fig[4, 1])
ax11 = Axis(fig[4, 2])
# ax12 = Axis(fig[4, 3])

# ax13 = Axis(fig[5, 1])
# ax14 = Axis(fig[5, 2])
# ax15 = Axis(fig[5, 3])

# ax16 = Axis(fig[6, 1])
# ax17 = Axis(fig[6, 2])
# ax18 = Axis(fig[6, 3])
ax = [
    ax1,
    ax2,
    # ax3,
    ax4,
    ax5,
    # ax6,
    ax7,
    ax8,
    # ax9,
    ax10,
    ax11,
    # ax12,
    # ax13,
    # ax14,
    # ax15,
    # ax16,
    # ax17,
    # ax18,
]
files = [
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.25_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.25_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.25_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.125_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.125_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.125_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.0625_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F-0.125_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.125_M1_d60_ΩTnothing_τ20.jld2",
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.125_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F-0.25_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.125_F-0.25_M1_d60_ΩTnothing_τ20.jld2",
    # "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.0625_F-0.25_M1_d60_ΩTnothing_τ20.jld2",
]

for ii = 1:length(ax)
    mkFigure(ax[ii], files[ii])

end

fig
save("test.pdf", fig)



## Single Trajectory Analysis
data = load_object(
    "data/Non_Thermal/Single_x0[5.5]_v0[0.5]_MemInf_s0.25_F0.03125_M1_d60_ΩTnothing_τ20.jld2",
)

step = 50


# l2 = lines!(ax, (1:length(Rs)) .* δ, [x[2] for x in res.Rs] .- 0, color = my_red)
# l2 = lines!(ax, (1:length(Rs)) .* δ, [x[3] for x in res.Rs] .- 0, color = my_green)
fig

# data.rs[1]



m = data.m                      # Mass of the chain atoms
k = data.k                      # Spring force constant
K = data.K                      # Confining potential force constant
Ωmin = sqrt(K / m)              # Smallest chain frequency
Ωmax = sqrt(4 * k / m + K / m)  # Largest chain frequency

M = data.M                      # Mass of the mobile particles
K_M = data.K_M                  # Trap force constant
ΩM = √(K_M / M)                 # Trap frequency
t_M = 2 * π / ΩM                # Period of the trapped mass
Rs = data.Rs |> vec
rs = data.rs
ts = data.ts ./ t_M

let
    fig = Figure(resolution = (1200, 600), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"t/t_M", ylabel = L"R")
    ax3 = Axis(fig[1, 2], xlabel = L"t/t_M", ylabel = L"R")
    ax5 = Axis(fig[1, 3], xlabel = L"t/t_M", ylabel = L"R")

    ax2 = Axis(fig[2, 1], xlabel = L"t/t_M", ylabel = L"r")
    ax4 = Axis(fig[2, 2], xlabel = L"t/t_M", ylabel = L"r")
    ax6 = Axis(fig[2, 3], xlabel = L"t/t_M", ylabel = L"r")

    ax_x_lims = [(0, 120), (0, 120), (88, 92), (88, 92), (116, 117), (116, 117)]
    ax_y_lims =
        [(-11, 11), (-0.15, 0.15), (-5.5, 5.5), (-0.15, 0.15), (-0.04, 0.04), (-0.04, 0.04)]
    ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    txt = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

    for n = 1:length(ax)
        xlims!(ax[n], ax_x_lims[n])
        ylims!(ax[n], ax_y_lims[n])
        text!(
            ax[n],
            txt[n],
            position = (
                ax_x_lims[n][1] + 5 / 6 * (ax_x_lims[n][2] - ax_x_lims[n][1]),
                ax_y_lims[n][2] / 2,
            ),
            textsize = 16,
            font = "CMU Serif",
        )
    end
    idx = findall(x -> (x <= ax_x_lims[1][2]) && (x >= ax_x_lims[1][1]), ts)
    lines!(ax[1], ts[idx], Rs[idx], color = my_red, linewidth = 1)
    lines!(ax[2], ts[idx], rs[idx], color = my_red, linewidth = 1)

    idx = findall(x -> (x <= ax_x_lims[3][2]) && (x >= ax_x_lims[3][1]), ts)
    lines!(ax[3], ts[idx], Rs[idx], color = my_red, linewidth = 1)
    lines!(ax[4], ts[idx], rs[idx], color = my_red, linewidth = 1)


    idx = findall(x -> (x <= ax_x_lims[5][2]) && (x >= ax_x_lims[5][1]), ts)
    lines!(ax[5], ts[idx], Rs[idx], color = my_red, linewidth = 1)
    lines!(ax[6], ts[idx], rs[idx], color = my_red, linewidth = 1)

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    save("General_Example.pdf", fig)
end

