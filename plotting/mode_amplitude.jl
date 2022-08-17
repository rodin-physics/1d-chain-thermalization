include("../src/main.jl")

data = load_object("data/non_thermal/Single_sigma0[220]_sigmadot0[40]_MemInf_lambda4_Phi-8_mu1_d60_bias0.1_omegaTnothing_tau200.jld2")
ρs = data.ρs
displacement = ρs .- data.α * (1 : 300)
function FT_Harm(k, disp)
    sum(exp.(1im * k .* (1:length(disp))) .* disp)
end
ks = range(0, 2 * pi, length = 100)
res = [FT_Harm(k, displacement[:, t]) for k in ks, t in 1:10:119999]
# res = [FT_Harm(k, displacement[:, t]) for k in ks, t in 1:size(displacement)[2]]
res = abs.(res)
heatmap(collect(1:20:119999), ks, res', colormap = :amp, title = L"\lambda = 4.0, \Phi = 4.0, speed = 40, bias = 0.2")
