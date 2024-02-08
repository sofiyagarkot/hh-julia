include("utils.jl")
using ForwardDiff

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# Colors
ioup_prior_color = "purple2"
ioup_transformed_prior_color = "green"
iwp_prior_color = "black"
iwp_transformed_prior_color = "olivedrab1"

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
                solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
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
        label = "IWP(3)",
        ylabel = "L2 loss", 
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
        label = "IOUP(3, rate)",
        ylabel = "L2 loss", 
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


# defining inverse and forward transforms
identity =  x -> 1*x
deriv_logit = x -> ForwardDiff.derivative(logit, x)
deriv_identity = x -> ForwardDiff.derivative(identity, x)

forward_transforms = [identity, logit, logit, logit]
deriv_forward_transforms = [deriv_identity, deriv_logit, deriv_logit, deriv_logit]
inverse_transforms = [identity, sigmoid, sigmoid, sigmoid]

# defining the problem
ode_func, u0, tspan, params = define_problem_and_parameters()

# transformed initial value problem
u0s_transformed = []
for i in 1:length(forward_transforms)
    forward_transform = forward_transforms[i]
    u0_transformed = forward_transform(u0[i])
    if isa(u0_transformed,  Taylor1{Float64})
        u0_transformed = inv_u(0)
    end
    push!(u0s_transformed, u0_transformed)
end
u0s_transformed = float.(u0s_transformed)

prob_transformed = ODEProblem(transformed_ode, u0s_transformed, evaltspan, params)

all_errors_t = []
for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        try
                solution = solve(prob_transformed, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        catch e
                push!(all_errors_t, NaN)
                continue
        else
                solution = solve(prob_transformed, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)

                reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution.t)
                
                # error = l2_transformed(hcat(solution.u...), reference(solution.t), inverse_transforms)
                error = l2_transformed(solution_transformed, reference, inverse_transforms)
                push!(all_errors_t, error)
        end
end


       
plot!( 
        p,
        rates, 
        all_errors_t,
        label = "Transformed IOUP(3, rate)",
        ylabel = "L2 loss", 
        guidefont=font(11, "Computer Modern"),
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
        label = "transformed IWP(3)",
        ylabel = "L2 loss", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        framestyle=:axes,
        fg_legend = :transparent
        )

display(p)

path = "./visuals/IOUP-rate-and-transformed.png"
savefig(p, path)


