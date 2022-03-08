# fig = Figure()
# ax = Axis(fig[1, 1])
# l1 = lines!(ax, res.ts, [x[1] for x in res.Rs])
# # l2 = lines!(ax, res.ts, [x[2] for x in res.Rs] , color = my_red)
# # l3 = lines!(ax, res.ts, [x[3] for x in res.Rs] , color = my_green)
# # l3 = lines!(ax, res.ts, [x[4] for x in res.Rs] , color = my_orange)
# for ii = 1:nChain
#     lines!(ax, res.ts, [x[ii] for x in res.rs], color = colorant"rgba(0, 0, 0, 0.35)")
# end
# # l2 = lines!(ax, (1:length(Rs)) .* δ, [x[2] for x in res.Rs] .- 0, color = my_red)
# # l2 = lines!(ax, (1:length(Rs)) .* δ, [x[3] for x in res.Rs] .- 0, color = my_green)
# fig
