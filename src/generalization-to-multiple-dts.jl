using RecursiveArrayTools
include("utils.jl")
using  LaTeXStrings

# Colors
ioup_prior_color = "purple4"
# ioup_transformed_prior_color = "green3"
iwp_prior_color = "black"
iwp_transformed_prior_color = "darkgreen"

dts = 10.0 .^ range(-1, -4, length=11)[1:end-1]
evaltspan = (0.0, 50.0)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()

rate_optimal_ioup = 12.9
# rate_optimal_ioup_transformed = 6.1

# the non-transformed priors

prob = define_problem()

p = plot(legendfont = font(7), size = (500, 400))

# IWP(3) 
prior = IWP(3)
algorithm = EK1(prior=prior, smooth=DENSE, diffusionmodel=DM)
errors_iwp = []
for dt in dts
    try

        solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)

        reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
        error = l2(solution.u, reference(solution.t))
        push!(errors_iwp, error)

    catch e
        
        push!(errors_iwp, NaN)
        continue

    else            
    end

end

plot!(p, dts, errors_iwp, label="IWP(3)", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_prior_color,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_prior_color,
    fg_legend = :transparent)
    


# IOUP(3) 

errors_ioup = []

for i in 1:length(dts)
    dt = dts[i]
    r = rate_optimal_ioup
    print("dt", dt, "\n")
    # rates = range(-300, 300, length=100)
    # r = optimize_rate(rates, prob, dt)
    # push!(rates_optimal, r)
    # print("rates_optimal", rates_optimal, "\n")
    algorithm = EK1(prior=IOUP(3,  r, update_rate_parameter=false), smooth=DENSE, diffusionmodel=DM)
    try
        solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)

        reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
        error = l2(solution.u, reference(solution.t))
        push!(errors_ioup, error)
    catch e
        push!(errors_ioup, NaN)
        continue
    
    end

end
# display(p)
plot!(p, dts[.!isnan.(errors_ioup)], 
    float.(errors_ioup[.!isnan.(errors_ioup)]), 
    linestyle=:dash, 
    label="IOUP(3, $rate_optimal_ioup)", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=ioup_prior_color,
    fg_legend = :transparent)

scatter!(p, dts[.!isnan.(errors_ioup)], 
    float.(errors_ioup[.!isnan.(errors_ioup)]), label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=ioup_prior_color,
    fg_legend = :transparent)

# the transformed priors
# IWP(3)

# function optimize_rate_transformed(rates, problem, problem_transformed, inverse_transforms, dt)
#     errors = []
#     for rate in rates
#         algorithm = EK1(prior=IOUP(3, rate, update_rate_parameter=false), smooth=DENSE, diffusionmodel=DM)
#         try
#             solution = solve(problem_transformed, algorithm, dt=dt, adaptive=false, dense = false)
#             reference = solve(problem, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution.t)
#             error = l2_transformed(solution, reference, inverse_transforms)
#             push!(errors, error)
#         catch e
#             print("error: ", e, "\n")
#             push!(errors, NaN)
#             continue
#         end
#     end
#     if all(isnan, errors)
#         return NaN
#     end
#     # min_error = minimum(errors[.!isnan.(errors)])
#     min_index = argmin(errors[.!isnan.(errors)])
#     p2 = plot(
#         rates[.!isnan.(errors)], 
#         float.(errors[.!isnan.(errors)]), 
#         label="IOUP(3), dt = $dt", 
#         framestyle=:axes,
#         # xlabel="rate", https://github.com/JuliaPlots/Plots.jl/issues/4816
#         # ylabel="L2 error", ERROR: AssertionError: total_plotarea_horizontal > 0mm
#         color=ioup_transformed_prior_color,
#         fg_legend = :transparent)

#     display(p2)

#     savefig(p2, "./visuals/link-funcs/rate_per_dt/dt-$dt.png")

#     return rates[min_index]
# end 

# transformed IWP(3)

prior = IWP(3)
prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()

errors_iwp_t = []
for dt in dts
    try
        solution_transformed = solve(prob_transformed,  EK1(prior=prior, smooth = SMOOTH, diffusionmodel=FixedDiffusion()), adaptive = false, dense=DENSE, dt=dt)
        reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
        error = l2_transformed(solution_transformed, reference, inverse_transforms)
        push!(errors_iwp_t, error)    
    catch e
        push!(errors_iwp_t, NaN)
        continue
    end
end



plot!(p, dts[.!isnan.(errors_iwp_t)], 
    errors_iwp_t[.!isnan.(errors_iwp_t)], 
    label="Transformed IWP(3)", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp_t, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color,
    fg_legend = :transparent)


# transformed IOUP(3)

# errors_ioup_t = []
# rates_optimal_transformed = []
# for dt in dts
#     print("dt", dt, "\n")
#     rates = range(-300, 300, length=100)
#     r = optimize_rate_transformed(rates, prob, prob_transformed, inverse_transforms, dt)
    
#     if isnan(r)
#         push!(errors_ioup_t, NaN)
#         continue
#     end

#     push!(rates_optimal_transformed, r)
#     print("rates_optimal_transformed", rates_optimal_transformed, "\n")
#     algorithm = EK1(prior=IOUP(3,  r, update_rate_parameter=false), smooth=DENSE, diffusionmodel=DM)
#     try
#         solution_transformed = solve(prob_transformed, algorithm,  dt=dt, adaptive=false, dense = false)
#         reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
#         error = l2_transformed(solution_transformed, reference, inverse_transforms)
#         push!(errors_ioup_t, error) 
#     catch e
#         push!(errors_ioup_t, NaN)
#         continue
#     end

# end

# plot!(p, dts[.!isnan.(errors_ioup_t)], 
#     errors_ioup_t[.!isnan.(errors_ioup_t)],
#     label="Transformed IOUP(3)", 
#     framestyle=:axes,
#     xaxis=:log10, yaxis=:log10,
#     color="darkred",
#     fg_legend = :transparent)

# scatter!(p, dts[.!isnan.(errors_ioup_t)],
#     errors_ioup_t[.!isnan.(errors_ioup_t)], label="", 
#     framestyle=:axes,
#     xaxis=:log10, yaxis=:log10,
#     color="darkred",
#     fg_legend = :transparent)
    


# savefig(p, "./visuals/multiple-dts-transformed-priors.png")
plot!(p, xaxis=:log10, yaxis=:log10, 
        legend=:bottomright, xlabel=L"$\Delta$t", 
        ylabel="L2 error",
        dpi=600)
p

savefig(p, "./visuals/presentation/multiple-dts-transformed-priors-num-3.png")