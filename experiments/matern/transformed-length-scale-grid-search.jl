include("../../src/utils.jl")

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

# Generate the log-spaced values for rate 
exponents = range(-6, stop=5, length=150)
lengthscales = 10 .^ exponents

matern_prior_transformed_color = "deepskyblue"
iwp_transformed_prior_color = "forestgreen"

algorithms = []
n_derivatives = 3
for l in lengthscales
        push!(algorithms, EK1(prior=Matern(n_derivatives, l), smooth=SMOOTH, diffusionmodel=DM))
end

all_errors = []
names_of_channels = ["V", "m", "h", "n"]

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

p = plot(layout=(1,1), legendfont = font(7), size = (500, 180))

xlabel = "length scale"

       
plot!( 
        p,
        lengthscales, 
        all_errors_t,
        label = "Transformed Matérn",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        legend=:topright,
        framestyle=:axes,
        fg_legend = :transparent,
        color= matern_prior_transformed_color
        )

prior = IWP(3)
solution_transformed = solve(prob_transformed,  EK1(prior=prior, smooth = SMOOTH, diffusionmodel=FixedDiffusion()), adaptive = false, dense=DENSE, tstops=tstops)
reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
error_iwp_t = l2_transformed(solution_transformed, reference, inverse_transforms)


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


plot!(p,  
        yaxis=:log10, 
        xlabel="length scale", 
        xaxis=:log10, 
        ylabel="tRMSE",
        dpi=600)
display(p)

savefig(p, "./visuals/matern/transformed-length-scale-grid-search.png")