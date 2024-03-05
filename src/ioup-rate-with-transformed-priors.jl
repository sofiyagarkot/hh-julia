include("utils.jl")
using ForwardDiff

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# Colors
ioup_prior_color = "purple4"
ioup_transformed_prior_color = "green3"
iwp_prior_color = "black"
iwp_transformed_prior_color = "darkgreen"

# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

# exponents = range(0, stop=40, length=30)
# rates = -10 .^ exponents
# rates = vcat(rates, -rates)
# rates = sort(rates)
rates = range(-30, stop=30, length=50)

# create the algorithms
algorithms = []
n_derivatives = 3
for l in rates
        push!(algorithms, EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false), smooth=SMOOTH, diffusionmodel=DM))
end

# calculate the errors
all_errors = []
names_of_channels = ["V", "m", "h", "n"]

for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        try
                solution = solve(prob, algorithm, dt=0.01, adaptive=ADAPTIVE, dense=DENSE)
        catch e
                push!(all_errors, NaN)
                continue
        else
                solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)

                reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution.t)
                
                error = l2(solution.u, reference(solution.t))
                push!(all_errors, error)
        end
end

# compare to baseline
solution_IWP = solve(prob, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM), tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_IWP.t)

error_IWP = l2(solution_IWP.u, reference(solution_IWP.t))


# plotting
p = plot(layout=(1,1), legendfont = font(7), size = (800, 180))

xlabel = "rate"

hline!( 
        p,
        [error_IWP],
        linestyle=:dash, linecolor=iwp_prior_color,
        label = "IWP",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        framestyle=:axes,
        fg_legend = :transparent
        )

plot!( 
        p,
        rates, 
        all_errors,
        label = "IOUP(rate)",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        # yticks=[],
        color = ioup_prior_color,
        framestyle=:axes,
        fg_legend = :transparent,
        )

prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()


all_errors_t = []
for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        try
                solution_transformed = solve(prob_transformed, algorithm, dt=0.01, adaptive=ADAPTIVE, dense=DENSE)
                reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
                
                # error = l2_transformed(hcat(solution.u...), reference(solution.t), inverse_transforms)
                error = l2_transformed(solution_transformed, reference, inverse_transforms)
                push!(all_errors_t, error)
        catch e
                push!(all_errors_t, NaN)
                continue
        end
end


       
plot!( 
        p,
        rates, 
        all_errors_t,
        label = "Transformed IOUP(rate)",
        ylabel = "tRMSE", 
        # guidefont=font(11, "Computer Modern"),
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        # yticks=[],
        framestyle=:axes,
        fg_legend = :transparent,
        color= ioup_transformed_prior_color
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
        label = "Transformed IWP(3)",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        framestyle=:axes,
        fg_legend = :transparent
        )


plot!(p,  yaxis=:log10, 
        # legend=:bottomright,
         xlabel="rate", 
        ylabel="L2 error",
        dpi=600)
display(p)

# path = "./visuals/IOUP-rate-and-transformed.png"
path = "./visuals/essay/IOUP-rate-and-transformed.png"
savefig(p, path)


