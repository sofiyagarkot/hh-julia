using RecursiveArrayTools
include("../../src/utils.jl")
using  LaTeXStrings

# Colors
iwp_prior_color = "black"
iwp_transformed_prior_color = "darkorange"
iwp_transformed_prior_color_V = "purple3"

iwp_transformed_prior_color_V_slope_2 = "springgreen"
iwp_transformed_prior_color_V_slope_3 = "limegreen"

dts = 10.0 .^ range(-1, -4, length=11)[1:end-1]
evaltspan = (0.0, 50.0)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()

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

plot!(p, dts, errors_iwp, label="IWP", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_prior_color,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_prior_color,
    fg_legend = :transparent)
    

#  transformed priors for channels
# IWP(3)

prior = IWP(3)
prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings()

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
    label="Transformed IWP\nOnly channels", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp_t, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color,
    fg_legend = :transparent)

# transformed priors for all the channels 
# IWP(3)

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



plot!(p,dts[.!isnan.(errors_iwp_t)], 
    errors_iwp_t[.!isnan.(errors_iwp_t)], 
    label="Transformed IWP\nV-scale = [-110, +60]", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color_V,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp_t, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color_V,
    fg_legend = :transparent)

# transformed priors for all the channels  with different slope 


function sigmoid_V(x, slope = 0.05, x_offset = 0.0, y_offset = -82.0, scale=115)
    return scale / (1 + exp(-slope * (x - x_offset))) + y_offset
end

function inverse_sigmoid_V(y, slope = 0.05, x_offset = 0.0, y_offset = -82.0, scale=115)
    return (-1/slope)*log( (scale)/(y - y_offset) - 1) + x_offset
end

# defining inverse and forward transforms
deriv_inverse_sigmoid_V =  x -> ForwardDiff.derivative(inverse_sigmoid_V, x)
deriv_logit = x -> ForwardDiff.derivative(logit, x)

forward_transforms = [inverse_sigmoid_V, logit, logit, logit]
deriv_forward_transforms = [deriv_inverse_sigmoid_V, deriv_logit, deriv_logit, deriv_logit]
inverse_transforms = [sigmoid_V, sigmoid, sigmoid, sigmoid]

# IWP(3)

prior = IWP(3)

prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V(forward_transforms, 
                                                                                                deriv_forward_transforms, inverse_transforms)           

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



plot!(p,dts[.!isnan.(errors_iwp_t)], 
    errors_iwp_t[.!isnan.(errors_iwp_t)], 
    label="Transformed IWP\nV-scale = [-82, +33]", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color_V_slope_2,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp_t, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_transformed_prior_color_V_slope_2,
    fg_legend = :transparent)


plot!(p, xaxis=:log10, yaxis=:log10, 
        dpi=600,
        legend=:bottomright, xlabel=L"$\Delta$t", 
        ylabel="tRMSE")
p

savefig(p, "./visuals/iwp/transformed-voltage-different-bounds.png")