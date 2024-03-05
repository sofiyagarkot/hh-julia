using ProbNumDiffEq, ComponentArrays, UnPack
using Plots, LaTeXStrings
using OrdinaryDiffEq
using Statistics, LinearAlgebra
using Plots.PlotMeasures
using Colors
using DiffEqDevTools
using Distributions
using ForwardDiff
using TaylorSeries
using RecursiveArrayTools

Plots.theme(
    :dao;
    markerstrokewidth=1.2,
    legend=:left,
)

function l2(sol, analytic)
    return sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
                                                      sol - analytic)))
end

function l2_transformed(sol_transformed, analytic, inferse_funcs)
    errors = []
    for i in 1:length(inferse_funcs)
        inferse_func = inferse_funcs[i]
        sol_channel = inferse_func.(sol_transformed[i,:])
        analytic_channel = analytic[i, :]
        error = sol_channel - analytic_channel
        push!(errors, error)
    end

    return sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2, vcat(errors...))))
end

function define_problem(duration = 50.0)
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

    function f(du, u, params, t)
        V, m, n, h = u

        @unpack gNa, gK = params

        I_inj = 500*1e-6

        du[2] = dmdt = (αm(V, VT) * (1 - m) - βm(V, VT) * m)
        du[3] = dndt = (αn(V, VT) * (1 - n) - βn(V, VT) * n)
        du[4] = dhdt = (αh(V, VT) * (1 - h) - βh(V, VT) * h)

        INa = gNa * m^3 * h * (V - ENa) * area
        IK = gK * n^4 * (V - EK) * area
        Ileak = gleak * (V - Eleak) * area
        Cm = C * area
        du[1] = dVdt = -(Ileak + INa + IK - I_inj) / Cm
    end

    u0 = [V0, m_inf(V0, VT), n_inf(V0, VT), h_inf(V0, VT)]
    tspan = (0.0, duration)
    p = ComponentVector(gNa=20.0, gK=15.0)

    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

function sigmoid(x, slope=1.0, x_offset = 0.0, y_offset= 0.0)
    return 1 / (1 + exp(-slope * (x - x_offset))) + y_offset
end

function logit(y, slope=1.0, x_offset=0.0, y_offset=0.0)
    """
    Inverse sigmoid function
    """
    return log((y - y_offset) / (1 - y + y_offset)) / slope + x_offset
end

function sigmoid_V(x, slope = 0.05, x_offset = 0.0, y_offset = -110.0, scale=170)
    return scale / (1 + exp(-slope * (x - x_offset))) + y_offset
end

function inverse_sigmoid_V(y, slope = 0.05, x_offset = 0.0, y_offset = -110.0, scale=170)
    return (-1/slope)*log( (scale)/(y - y_offset) - 1) + x_offset
end

# function inverse_sigmoid_V(y, slope = 100.0, x_offset = 0.0, y_offset = -110.0, scale=170)
#     return (1/slope)*log( -(y-y_offset)/(y - (scale + y_offset) )) + x_offset
# end

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

function generate_transformed_settings(forward_transforms= [], deriv_forward_transforms= [], inverse_transforms= [])
    
    if isempty(forward_transforms)
        
        # defining inverse and forward transforms
        identity =  x -> 1*x
        deriv_logit = x -> ForwardDiff.derivative(logit, x)
        deriv_identity = x -> ForwardDiff.derivative(identity, x)

        forward_transforms = [identity, logit, logit, logit]
        deriv_forward_transforms = [deriv_identity, deriv_logit, deriv_logit, deriv_logit]
        inverse_transforms = [identity, sigmoid, sigmoid, sigmoid]

    end
    
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

    return prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func
end 



function generate_transformed_settings_V(forward_transforms= [], deriv_forward_transforms= [], inverse_transforms= [])


    if isempty(forward_transforms)
        
        
        # defining inverse and forward transforms
        deriv_inverse_sigmoid_V =  x -> ForwardDiff.derivative(inverse_sigmoid_V, x)
        deriv_logit = x -> ForwardDiff.derivative(logit, x)

        forward_transforms = [inverse_sigmoid_V, logit, logit, logit]
        deriv_forward_transforms = [deriv_inverse_sigmoid_V, deriv_logit, deriv_logit, deriv_logit]
        inverse_transforms = [sigmoid_V, sigmoid, sigmoid, sigmoid]

    end
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

    return prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func
end 





# function l2_in_time(sol, analytic)
#     return  map(x ->sqrt.( x .^ 2),  sol - analytic)
# end


# function l2_transformed_in_time(sol, analytic)
#     errors = []
#     for i in 1:length(inferse_funcs)
#         inferse_func = inferse_funcs[i]
#         sol_channel = inferse_func.(sol_transformed[i,:])
#         analytic_channel = analytic[i, :]
#         error = map(x ->sqrt.( x .^ 2), sol_channel - analytic_channel)
#         push!(errors, error)
#     end

#     return  errors
# end