using DifferentialEquations 
using ProbNumDiffEq
using Plots
using ForwardDiff
using TaylorSeries
using Statistics,LinearAlgebra
using Plots.PlotMeasures

function logistic_ode(u_vector, p, t)
    u = u_vector[1]

    du = u*(1-u)
    [du]
end

u0 = [0.05]
tspan = (0.0, 10.0)
prob = ODEProblem(logistic_ode, u0, tspan)
solution = solve(prob,  Tsit5())
plot(solution)

solution = solve(prob,  EK1(prior=IWP(3), smooth = true, diffusionmodel=FixedDiffusion()))
plot(solution)

function transformed_logistic_ode(u_vector, p , t )
    u = u_vector[1]
    du = 1-exp(u)
    [du]
end
u0_transformed = log.(u0)
prob_transformed = ODEProblem(transformed_logistic_ode, u0_transformed, tspan)
tstops = range(tspan..., length=100)

solution = solve(prob_transformed,  EK1(prior=IWP(3), smooth = true, diffusionmodel=FixedDiffusion()), dense=true, tstops=tstops)
nsamples = 100
samples = ProbNumDiffEq.sample(solution, nsamples)


p = plot(layout = (1, 1), legendfont = font(7), size = (800, 300))
for i in 1:nsamples
    plot!(p, solution.t, exp.(samples[:,:,i]), label="")
end


function transform_sigmoid(x, slope=1.0, x_offset = 0.0, y_offset= 0.0)
    return 1 / (1 + exp(-slope * (x - x_offset))) + y_offset
end

function inverse_logit(y, slope=1.0, x_offset=0.0, y_offset=0.0)
    """
    Inverse sigmoid function
    """
    return log((y - y_offset) / (1 - y + y_offset)) / slope + x_offset
end

function plot_transformed(
    forward_transform, deriv_forward_transform, inverse_transform, ode_func,
    u0, tspan; title="", prior=IWP(3), nsamples=200, tspan_extrapolate = (0.0, 11.0)
    )

    tstops=range(tspan..., length=Int(tspan[2]*100))
    u0_transformed = forward_transform.(u0)
    function transformed_ode(u_vector, p, t)
        dus =[]
        for i in 1:length(u_vector)
            u = u_vector[i]
            inv_u = inverse_transform(u)
            if isa(inv_u,  Taylor1{Float64})
                inv_u = inv_u(0)
            end
            deriv = deriv_forward_transform(inv_u)
            # deriv = deriv_forward_transform(u)
            
            du = deriv*ode_func(inv_u, p, t)[i]
            push!(dus, du)
        end
        return dus
    end
    
    
    prob_transformed = ODEProblem(transformed_ode, u0_transformed, tspan)
    solution = solve(prob_transformed,  EK1(prior=prior, smooth = true, diffusionmodel=FixedDiffusion()), dense=true, tstops=tstops)
    
    # OPTION 1 
    # from Nathanael's banner
    ts = range(tspan_extrapolate..., length=Int(tspan_extrapolate[2]*100))

    interp(sol, ts) = StructArray([
        ProbNumDiffEq.interpolate(
            t,
            sol.t,
            sol.x_filt,
            sol.x_smooth,
            sol.diffusions,
            sol.cache;
            smoothed=sol.alg.smooth,
        )
        for t in ts
    ])

    derivative = 0
    p1 = plot(layout = (length(u0_transformed), 1), legendfont = font(7), size = (1000, 600))
    for ch in 1:length(u0_transformed)
        xs = interp(solution, ts)

        vecvec2mat(vv) = hcat(vv...)'
        
        if derivative == 0
            H = solution.cache.E0
        elseif derivative == 1
            H = solution.cache.E1
        elseif derivative == 2
            H = solution.cache.E2
        elseif derivative == 3
            H = solution.cache.E3
        end
            
        title1 = "Mean and variance"
        
        m = vecvec2mat(xs.μ) * H'
        std__ = sqrt.(vecvec2mat(diag.(xs.Σ))) * H'
        plot!(p1[ch], ts, inverse_transform.(m[:, ch]), 
                    fillalpha = 0.5, xlabel = "t", ylabel = "u(t)", 
                    title=title1, legend=:topleft, label="", 
                    left_margin = [20mm 0mm], right_margin = [20mm 0mm])
        plot!(p1[ch], ts, inverse_transform.(m[:, ch] - 3*std__[:, ch]),  label="mean - 3 std",linestyle=:dash, title=title1)
        plot!(p1[ch], ts, inverse_transform.(m[:, ch] + 3*std__[:, ch]),  label="mean + 3 std", linestyle=:dash, title=title1)
                
    end
    display(p1)
    
    samples = ProbNumDiffEq.sample(solution, nsamples)

    p = plot(layout = (1, 1), legendfont = font(7), size = (800, 300))
    
    m = mean(samples, dims=3)
    std_ = std(samples, dims=3)
    scale = 3
    lower_bound = m - scale*std_
    upper_bound = m + scale*std_

    plot!(p, solution.t, inverse_transform.(m[:,:,1]),  label="", title=title)
    plot!(p, solution.t, inverse_transform.(lower_bound[:, :, 1]),  label="mean - $scale std",linestyle=:dash, title=title)
    plot!(p, solution.t, inverse_transform.(upper_bound[:, :, 1]),  label="mean + $scale std", linestyle=:dash, title=title)

    display(p)
    return p
end

forward_transform = log
inverse_transform = exp
deriv_forward_transform = x -> ForwardDiff.derivative(forward_transform, x)
tspan = (0.0, 1.0)
p = plot_transformed(
    forward_transform, deriv_forward_transform, inverse_transform,
    logistic_ode, u0, tspan ; title="Forward = exp, backward=log",
    tspan_extrapolate = (0.0, 5.5)
)

forward_transform = inverse_logit 
inverse_transform = transform_sigmoid
deriv_forward_transform = x -> ForwardDiff.derivative(forward_transform, x)
tspan = (0.0, 1.0)
p = plot_transformed(
    forward_transform, deriv_forward_transform, inverse_transform,
    logistic_ode, u0, tspan ; title="Forward = sigmoid, backward=inverse_logit",
    tspan_extrapolate = (0.0, 6.5)
)

forward_transform =  x -> 1*x
inverse_transform = x -> 1*x
deriv_forward_transform = x -> ForwardDiff.derivative(forward_transform, x)
tspan = (0.0, 1.0)
p = plot_transformed(
    transform_sigmoid, deriv_forward_transform, inverse_transform,
    logistic_ode, u0, tspan ; title="Forward = I, backward=I", tspan_extrapolate = (0.0, 6.5)
)