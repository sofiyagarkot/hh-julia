using DifferentialEquations
using Plots
using UnPack
using ComponentArrays
using ProbNumDiffEq
using SciMLBase
using LinearAlgebra
using Optimization
using Optim
using OptimizationOptimJL
using OrdinaryDiffEq
using Statistics



αm(V, VT) = -0.32 * (V - VT - 13) / (exp(-(V - VT - 13) / 4) - 1)
βm(V, VT) = 0.28 * (V - VT - 40) / (exp((V - VT - 40) / 5) - 1)

αn(V, VT) = -0.032 * (V - VT - 15) / (exp(-(V - VT - 15) / 5) - 1)
βn(V, VT) = 0.5 * exp(-(V - VT - 10) / 40)

αh(V, VT) = 0.128 * exp(-(V - VT - 17) / 18)
βh(V, VT) = 4 / (1 + exp(-(V - VT - 40) / 5))

# Would be the solution when dm/dt = 0:
m_inf(V, VT) = 1 / (1 + βm(V, VT) / αm(V, VT))
n_inf(V, VT) = 1 / (1 + βn(V, VT) / αn(V, VT))
h_inf(V, VT) = 1 / (1 + βh(V, VT) / αh(V, VT))

# Hodgkin-Huxley (Pospischil) with 2 free parameters
ENa = 53
EK = -107
area = 15e-5
C = 1
Eleak = -70
VT = -60
gleak = 0.1
V0 = -70

I(t::ProbNumDiffEq.Taylor1) = zero(t)
I(t::ProbNumDiffEq.TaylorN) = zero(t)

I(t::Float64) = (10 <= t <= 100) ? 500one(t) : zero(t)

function f(du, u, params, t)
    V, m, n, h = u

    @unpack gNa, gK = params

    I_inj = I(t) * 1e-6 # uA

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
tspan = (0.0, 110.0)
p = ComponentVector(gNa=20.0, gK=15.0)
prob = ODEProblem(f, u0, tspan, p)

# Matern prior
sol_matern = solve(prob, RosenbrockExpEK(prior=Matern(3, 1), smooth=true), dt=0.01, adaptive=false)
sol_ioup = solve(prob, RosenbrockExpEK(prior=IOUP(3, update_rate_parameter=true), smooth=true), dt=0.01, adaptive=false)
sol_iwp = solve(prob, RosenbrockExpEK(prior=IWP(3), smooth=true), dt=0.01, adaptive=false)
reference_matern = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_matern.t)


# _setups = [
#   "EK0(2)" => Dict(:alg=>EK0(order=2, smooth=DENSE))
#   "EK0(3)" => Dict(:alg=>EK0(order=3, smooth=DENSE))
#   "EK1(2)" => Dict(:alg=>EK1(order=2, smooth=DENSE))
#   "EK1(3)" => Dict(:alg=>EK1(order=3, smooth=DENSE))
#   "EK1(5)" => Dict(:alg=>EK1(order=5, smooth=DENSE))
#   "RosenbrockExpEK1(3)" => Dict(:alg=>RosenbrockExpEK(order=3, smooth=DENSE))
#   "RosenbrockExpEK1(5)" => Dict(:alg=>RosenbrockExpEK(order=5, smooth=DENSE))
# ]

errors_matern = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_matern.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_matern.pu))[ch, :]
     errors_ = ((abs.(sol_ .- ref_) )./ ref_) 
     push!(errors_matern, errors_)
end


# errors_matern = reduce(hcat, (mean.(sol_matern.pu) .- reference_matern.u)./ reference_matern.u)' 
# errors_matern = reduce(hcat, mean.(sol_matern.pu) .- reference_matern.u)' 
# errors_matern = reduce(hcat, (mean.(sol_matern.pu) .- reference_matern.u) ./ reference_matern.u)'

reference_iwp = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_iwp.t)

errors_iwp = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_iwp.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_iwp.pu))[ch, :]
     errors_ = ((abs.(sol_ .- ref_) )./ ref_) 
     push!(errors_iwp, errors_)
end

# errors_iwp = reduce(hcat, mean.(sol_iwp.pu) .- reference_iwp.u)' ./ reference_iwp.u

reference_ioup = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_ioup.t)
errors_ioup = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_ioup.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_ioup.pu))[ch, :]
     errors_ = ((abs.(sol_ .- ref_) )./ ref_) 
     push!(errors_ioup, errors_)
end

# errors_ioup = reduce(hcat, mean.(sol_ioup.pu) .- reference_ioup.u)' ./reference_ioup.u

using Plots.PlotMeasures
p1 = plot(sol_matern.t, errors_matern[1], label="Matern", xlabel="t", ylabel="error (V), %", title="RosenbrockExpEK", left_margin = [20mm 0mm], upper_margin = 40px, right_margin = [20mm 0mm])
p2 = plot(sol_matern.t, errors_matern[2], label="Matern", xlabel="t", ylabel="error (m), %",  left_margin = [20mm 0mm], right_margin = [20mm 0mm])
p3 = plot(sol_matern.t, errors_matern[3], label="Matern", xlabel="t", ylabel="error (n), %",  left_margin = [20mm 0mm], right_margin = [20mm 0mm])
p4 = plot(sol_matern.t, errors_matern[4], label="Matern", xlabel="t", ylabel="error (h), %",  left_margin = [20mm 0mm], bottom_margin = 40px, right_margin = [20mm 0mm])

plot!(p1, sol_iwp.t, errors_iwp[1],label="IWP", xlabel="t")
plot!(p1, sol_ioup.t, errors_ioup[1],label="IOUP", xlabel="t")

plot!(p2, sol_iwp.t, errors_iwp[2],label="IWP", xlabel="t")
plot!(p2, sol_ioup.t, errors_ioup[2],label="IOUP", xlabel="t")

plot!(p3, sol_iwp.t, errors_iwp[3],label="IWP", xlabel="t")
plot!(p3, sol_ioup.t, errors_ioup[3],label="IOUP", xlabel="t")


plot!(p4, sol_iwp.t, errors_iwp[4],label="IWP", xlabel="t")
plot!(p4, sol_ioup.t, errors_ioup[4],label="IOUP", xlabel="t")

l = @layout [a ; b ; c ; d]
plot(p1, p2, p3, p4, layout = l)
plot!(size=(1000,600))

# plot(errors_iwp,
#      legend=false,
#      layout=(4,1),
#      title=["Hodgkin-Huxley Errors (IWP prior, non-adaptive steps)" "" "" ""],
#      ylabel=["V(t)" "m(t)" "n(t)" "h(t)"],
#      xlabel=["" "" "" "t"],
#      size = (1000, 600),
#      color=[1 2 3 4],
#      xticks=:auto, yticks=:auto
# )


# plot(errors_ioup,
#      legend=false,
#      layout=(4,1),
#      title=["Hodgkin-Huxley Errors (IOUP prior, non-adaptive steps)" "" "" ""],
#      ylabel=["V(t)" "m(t)" "n(t)" "h(t)"],
#      xlabel=["" "" "" "t"],
#      size = (1000, 600),
#      color=[1 2 3 4],
#      xticks=:auto, yticks=:auto
# )


# plot(errors_matern,
#      legend=false,
#      layout=(4,1),
#      title=["Hodgkin-Huxley Errors (Matern prior, non-adaptive steps)" "" "" ""],
#      ylabel=["V(t)" "m(t)" "n(t)" "h(t)"],
#      xlabel=["" "" "" "t"],
#      size = (1000, 600),
#      color=[1 2 3 4],
#      xticks=:auto, yticks=:auto
# )