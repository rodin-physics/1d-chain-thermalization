include("../../src/main.jl")
## Parameters
α = 10
μ = 1
ωmax = 10
Φ0 = 5.0
λ = 1.0
ωT = 0.0

function speed_after_pass(init_speed)
    pred_mean = Δ_thermal_analytic(init_speed, Φ0, λ, ωmax, ωT)
    pred_var = Δ_thermal_variance(init_speed, Φ0, λ, ωmax, ωT)

    dist = Normal(pred_mean, √(pred_var))
    

end

function random_walk(init_speed, numP, nPasses)
    final_speeds = repeat([init_speed], numP)



end