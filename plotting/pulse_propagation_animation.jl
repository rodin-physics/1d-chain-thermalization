include("../src/main.jl")

s = load_object("precomputed/systems/System_K1_k20_m1_d60_l500.jld2")

Gs = s.G

@showprogress for jj = 0:4000
    fig = Figure(resolution = (1200, 1200), font = "CMU Serif", fontsize = 14)

    ax1 = Axis(fig[1, 1], xlabel = L"j", ylabel = L"\Delta r")
    sc = scatter!(ax1, 0:500, Gs[end-10*jj], markersize = 2)
    t_fast = round(10 * jj / 60, digits = 2)
    t_slow = round(10 * jj / 60 / 9, digits = 2)
    text!(
        ax1,
        L"t = %$t_fast \times t_{fast} = %$t_slow \times t_{slow}",
        position = (400, 0.1),
        textsize = 16,
        font = "CMU Serif",
    )
    ylims!(ax1, (-0.15, 0.15))
    save("anim/frame$(jj).png", fig)
end
fig
