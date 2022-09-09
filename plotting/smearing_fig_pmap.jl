using Distributed
proc_number = 25
if nprocs() < proc_number
    addprocs(proc_number - nprocs())
end
@everywhere using DelimitedFiles
@everywhere include("smearing_fig.jl")


Δ_bias_rep = @showprogress pmap(x -> Δ_transport_smeared(x, Φ0, λ, ωmax, α, μ), xs2)
Δ_bias_att = @showprogress pmap(x -> Δ_transport_smeared(x, -Φ0, λ, ωmax, α, μ), xs)

writedlm("rep_Deltabar.dat", Δ_bias_rep)
writedlm("att_Deltabar.dat", Δ_bias_att)
