include("../../src/utils.jl")
using ForwardDiff

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# Colors
ioup_prior_color = "purple4"
ioup_transformed_prior_color = "orange"
iwp_prior_color = "black"
iwp_transformed_prior_color = "forestgreen"

# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)
rates = range(-30, stop=30, length=50)

# create the algorithms
algorithms = []
n_derivatives = 3
for l in rates
        push!(algorithms, EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false), smooth=SMOOTH, diffusionmodel=DM))
end

prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()


all_errors_t = []
for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        try
                solution_transformed = solve(prob_transformed, algorithm, dt=0.01, adaptive=ADAPTIVE, dense=DENSE)
                reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
                error = l2_transformed(solution_transformed, reference, inverse_transforms)
                push!(all_errors_t, error)
        catch e
                push!(all_errors_t, NaN)
                continue
        end
end


prior = IWP(3)
solution_transformed = solve(prob_transformed,  EK1(prior=prior, smooth = SMOOTH, diffusionmodel=FixedDiffusion()), adaptive = false, dense=DENSE, tstops=tstops)
reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
error_iwp_t = l2_transformed(solution_transformed, reference, inverse_transforms)


# plotting
p = plot(layout=(1,1), legendfont = font(7), size = (500, 180))

xlabel = "rate"
       
plot!( 
        p,
        rates, 
        all_errors_t,
        label = "Transformed IOUP",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        legend=:topright,
        framestyle=:axes,
        fg_legend = :transparent,
        color= ioup_transformed_prior_color
        )

hline!( 
        p,
        [error_iwp_t],
        linestyle=:dash, 
        linecolor=iwp_transformed_prior_color,
        label = "Transformed IWP",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        left_margin = [7mm 0mm], 
        right_margin = [1mm 0mm],
        bottom_margin = [7mm 0mm],
        legend=:topright,
        framestyle=:axes,
        fg_legend = :transparent
        )


plot!(p,  yaxis=:log10, 
         xlabel="rate", 
        ylabel="tRMSE",
        dpi=600)
display(p)

path = "./visuals/ioup/transformed-rate-grid-search.png"
savefig(p, path)


