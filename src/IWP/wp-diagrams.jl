using LinearAlgebra, Statistics
using Plots.PlotMeasures
using ProbNumDiffEq
using ComponentArrays
using LinearAlgebra, Statistics
using UnPack
using Plots 
using StaticArrays, DiffEqDevTools, ParameterizedFunctions, Plots, SciMLBase, OrdinaryDiffEq

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

DENSE = SMOOTH = TO_SAVE = false

dts = 10.0 .^ range(-2, -6, length=10)[begin:end-1]
abstols = reltols = repeat([missing], length(dts))

ADAPTIVE = false

DM = FixedDiffusion()

_setups = [
    # "DP5" => Dict(:alg=>DP5(), :dts=>dts)
    # "Tsit5" => Dict(:alg=>Tsit5(), :dts=>dts)
    # "Vern7" => Dict(:alg=>Vern7(), :dts=>dts)
    "Exponential Euler" => Dict(:alg=>LawsonEuler(krylov = true, m = 50), :dts=>dts)
    "EK1(1)" => Dict(:alg=>EK1(order=1, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    "EK1(2)" => Dict(:alg=>EK1(order=2, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    "EK1(3)" => Dict(:alg=>EK1(order=3, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    "EK1(4)" => Dict(:alg=>EK1(order=4, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    ]
    
# print("setups", _setups)
labels = first.(_setups)
setups = last.(_setups)

prob = define_problem()
test_sol = solve(prob, Vern9(), abstol=1/10^14, reltol=1/10^14)
SAVE_EVERYSTEP = false


wp = WorkPrecisionSet(
        prob, abstols, reltols, setups;
        names = labels,
        appxsol = test_sol,
        dense = DENSE,
        save_everystep = SAVE_EVERYSTEP,
        maxiters = Int(1e7),
        numruns = 5,
        error_estimate = :l2,
        adaptive = ADAPTIVE,
    )
colors = [colorant"darkorange"]

c1 = colorant"lightblue"
c2 = colorant"darkblue"
colors2 = range(c1, stop=c2, length=4)
# push!(colors, colorant"darkorange")
for c in colors2
    push!(colors, c)
end
p = plot(wp, title="", color=colors', legend=:bottomright, legendfontsize=8)

p


savefig(p, "./visuals/baseline/fixed_diffusion_wp_EK1_IWP_ExpEuler.png")

p2 = plot(wp, x=:dts, title="", 
        fontsize=7, 
        legend=:bottomright, 
        color=colors',
        legendfontsize=8,
        )
p2

savefig(p2, "./visuals/baseline/fixed_diffusion_steps_number_wp_EK1_IWP_ExpEuler.png")


p3 = plot(wp, x=:dts,y=:l2, title="", 
        fontsize=7, 
        legend=:bottomright, 
        color=colors',
        legendfontsize=8,
        )
p3
savefig(p3, "./visuals/baseline/fixed_diffusion_steps_number_l2_error_wp_EK1_IWP_ExpEuler.png")
