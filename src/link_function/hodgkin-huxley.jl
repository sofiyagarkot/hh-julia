using ForwardDiff
using TaylorSeries
using LinearAlgebra, Statistics
using StaticArrays, DiffEqDevTools, ParameterizedFunctions, Plots, SciMLBase, OrdinaryDiffEq
using Plots.PlotMeasures
using ProbNumDiffEq
using ComponentArrays

function define_problem_and_parameters(duration = 50.0)
    αm(V, VT) = -0.32 * (V - VT - 13) / (exp(-(V - VT - 13) / 4) - 1)
    βm(V, VT) = 0.28 * (V - VT - 40) / (exp((V - VT - 40) / 5) - 1)

    αn(V, VT) = -0.032 * (V - VT - 15) / (exp(-(V - VT - 15) / 5) - 1)
    βn(V, VT) = 0.5 * exp(-(V - VT - 10) / 40)

    αh(V, VT) = 0.128 * exp(-(V - VT - 17) / 18)
    βh(V, VT) = 4 / (1 + exp(-(V - VT - 40) / 5))

    m_inf(V, VT) = 1 / (1 + βm(V, VT) / αm(V, VT))
    n_inf(V, VT) = 1 / (1 + βn(V, VT) / αn(V, VT))
    h_inf(V, VT) = 1 / (1 + βh(V, VT) / αh(V, VT))

    ENa = 53
    EK = -107
    area = 15e-5
    C = 1
    Eleak = -70
    VT = -60
    gleak = 0.1
    V0 = -70

    I(t) = (0 <= t <= duration) ? 500one(t) : zero(t)

    function f(u, params, t)
        
        V, m, n, h = u

        @unpack gNa, gK = params

        I_inj = 500*1e-6

        dmdt = (αm(V, VT) * (1 - m) - βm(V, VT) * m)
        dndt = (αn(V, VT) * (1 - n) - βn(V, VT) * n)
        dhdt = (αh(V, VT) * (1 - h) - βh(V, VT) * h)

        INa = gNa * m^3 * h * (V - ENa) * area
        IK = gK * n^4 * (V - EK) * area
        Ileak = gleak * (V - Eleak) * area
        Cm = C * area
        dVdt = -(Ileak + INa + IK - I_inj) / Cm

        [dVdt, dmdt, dndt, dhdt]
    end

    u0 = [V0, m_inf(V0, VT), n_inf(V0, VT), h_inf(V0, VT)]
    tspan = (0.0, duration)
    p = ComponentVector(gNa=20.0, gK=15.0)

    return f, u0, tspan, p
end

ode_func, u0, tspan, p = define_problem_and_parameters()

function sigmoid(x, slope=1.0, x_offset = 0.0, y_offset= 0.0)
    return 1 / (1 + exp(-slope * (x - x_offset))) + y_offset
end

function logit(y, slope=1.0, x_offset=0.0, y_offset=0.0)
    """
    Inverse sigmoid function
    """
    return log((y - y_offset) / (1 - y + y_offset)) / slope + x_offset
end

identity =  x -> 1*x

deriv_logit = x -> ForwardDiff.derivative(logit, x)
deriv_identity = x -> ForwardDiff.derivative(identity, x)

forward_transforms = [identity, logit, logit, logit]
deriv_forward_transforms = [deriv_identity, deriv_logit, deriv_logit, deriv_logit]
inverse_transforms = [identity, sigmoid, sigmoid, sigmoid]


function plot_transformed(
    forward_transforms, 
    deriv_forward_transforms, 
    inverse_transforms, 
    ode_func,
    p,
    u0, tspan; 
    titles=["","","", ""],
    xlabels=["","","", "t"],
    ylabels=["V","m","n", "h"],
    legends=[false, false, false, true],
    prior=IWP(3), 
    tspan_extrapolate = (0.0, 11.0)
    )

    tstops=range(tspan..., length=Int(tspan[2]*100))

    u0s_transformed = []
    for i in 1:length(forward_transforms)
        forward_transform = forward_transforms[i]
        u0_transformed = forward_transform(u0[i])
        if isa(u0_transformed,  Taylor1{Float64})
            u0_transformed = inv_u(0)
        end
        push!(u0s_transformed, u0_transformed)
    end

    function transformed_ode(u_vector, p, t)

        dus =[]

        inv_us = []
        derivs = []
        for i in 1:length(u_vector)   #  loop over all the ODE dimensions
            
            deriv_forward_transform = deriv_forward_transforms[i]
            inverse_transform  = inverse_transforms[i]

            u = u_vector[i]
            inv_u = inverse_transform(u)
            push!(inv_us, inv_u)

            if isa(inv_u,  Taylor1{Float64})
                inv_u = inv_u(0)
            end
            deriv = deriv_forward_transform(inv_u)
            push!(derivs, deriv)
        end

        dus = derivs.*ode_func(inv_us, p, t)
        return dus
    end
    
    u0s_transformed = float.(u0s_transformed)

    prob_transformed = ODEProblem(transformed_ode, u0s_transformed, tspan, p)
    solution = solve(prob_transformed,  EK1(prior=prior, smooth = true, diffusionmodel=FixedDiffusion()), dense=true, tstops=tstops)
    
    # OPTION 1 
    # from Nathanael's banner
    ts = range(tspan_extrapolate..., length=Int(round(tspan_extrapolate[2]*100)))

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
    p1 = plot(layout = (length(u0s_transformed), 1), legendfont = font(7), size = (1000, 600))
    for ch in 1:length(u0s_transformed)
        inverse_transform  = inverse_transforms[ch]
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
    
        
        m = vecvec2mat(xs.μ) * H'
        std__ = sqrt.(vecvec2mat(diag.(xs.Σ))) * H'
        plot!(p1[ch], ts, inverse_transform.(m[:, ch]), 
                    fillalpha = 0.5, 
                    xlabel = xlabels[ch], 
                    ylabel = ylabels[ch], 
                    title=titles[ch], 
                    legend=legends[ch], 
                    label="mean", 
                    left_margin = [20mm 0mm], 
                    right_margin = [20mm 0mm])
        plot!(p1[ch], ts, inverse_transform.(m[:, ch] - 3*std__[:, ch]),  
                    label="mean - 3 std",
                    linestyle=:dash, 
                    xlabel = xlabels[ch], 
                    ylabel = ylabels[ch], 
                    title=titles[ch], 
                    legend=legends[ch])
        plot!(p1[ch], ts, inverse_transform.(m[:, ch] + 3*std__[:, ch]),  label="mean + 3 std", linestyle=:dash, 
                    xlabel = xlabels[ch], 
                    ylabel = ylabels[ch], 
                    title=titles[ch], 
                    legend=legends[ch], )
                
    end
    display(p1)

    # OPTION 2
    # estimating mean and variance by sampling k samples from the solution
    
    # samples = ProbNumDiffEq.sample(solution, nsamples)

    # p = plot(layout = (1, 1), legendfont = font(7), size = (800, 300))
    
    # m = mean(samples, dims=3)
    # std_ = std(samples, dims=3)
    # scale = 3
    # lower_bound = m - scale*std_
    # upper_bound = m + scale*std_

    # plot!(p, solution.t, inverse_transform.(m[:,:,1]),  label="", title=title)
    # plot!(p, solution.t, inverse_transform.(lower_bound[:, :, 1]),  label="mean - $scale std",linestyle=:dash, title=title)
    # plot!(p, solution.t, inverse_transform.(upper_bound[:, :, 1]),  label="mean + $scale std", linestyle=:dash, title=title)

    # display(p)
    # return p
end

plot_transformed(
    forward_transforms,
    deriv_forward_transforms, 
    inverse_transforms,
    ode_func, p,
    u0, tspan ;
    titles=["Sigmoid on channels, and identity on Voltage", "", "", ""],
    tspan_extrapolate = (0.0, 52.3)
)