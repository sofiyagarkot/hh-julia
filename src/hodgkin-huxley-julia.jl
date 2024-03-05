module Hodgkin-Huxley-Julia

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

export define_problem

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