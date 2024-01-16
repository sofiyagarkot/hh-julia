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

sol_iwp3 = solve(prob, RosenbrockExpEK(prior=IWP(3), smooth=true), dt=0.01, adaptive=false)
sol_iwp5 = solve(prob, EK1(prior=IWP(3), smooth=true), dt=0.01, adaptive=false)
sol_matern3 = solve(prob, RosenbrockExpEK(prior=Matern(3, 1), smooth=true), dt=0.01, adaptive=false)
sol_matern5 = solve(prob, EK1(prior=Matern(3, 1), smooth=true), dt=0.01, adaptive=false)
sol_ioup3 = solve(prob, RosenbrockExpEK(prior=IOUP(3, update_rate_parameter=true), smooth=true), dt=0.01, adaptive=false)
sol_ioup5 = solve(prob, EK1(prior=IOUP(3, update_rate_parameter=true), smooth=true), dt=0.01, adaptive=false)

# sol_iwp3 = solve(prob, RosenbrockExpEK(prior=IWP(3), smooth=true), dt=0.01, adaptive=false)
# sol_iwp5 = solve(prob, EK1(prior=IWP(3), smooth=true), dt=0.01, adaptive=false)


reference_iwp3 = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_iwp3.t)
errors_iwp3 = []
stds_iwp3 = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_iwp3.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_iwp3.pu))[ch, :]
     errors_ = ((abs.(sol_ .- ref_) )./ ref_) .*100 
     sol_std = reduce(hcat, std.(sol_iwp3.pu))[ch, :]
     errors_ = (abs.((sol_ .- ref_) ./ ref_) + 3sol_std) .*100 
     push!(stds_iwp3, sol_std)
     push!(errors_iwp3, errors_)
end

reference_iwp5 = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_iwp5.t)
errors_iwp5 = []
stds_iwp5 = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_iwp5.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_iwp5.pu))[ch, :]
     sol_std = reduce(hcat, std.(sol_iwp5.pu))[ch, :]
     errors_ = (abs.((sol_ .- ref_) ./ ref_) + 3sol_std) .*100 
     push!(stds_iwp5, sol_std)
     push!(errors_iwp5, errors_)
end


reference_matern3 = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_matern3.t)
errors_matern3 = []
stds_matern3 = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_matern3.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_matern3.pu))[ch, :]
     errors_ = ((abs.(sol_ .- ref_) )./ ref_) 
     sol_std = reduce(hcat, std.(sol_matern3.pu))[ch, :]
     errors_ = (abs.((sol_ .- ref_) ./ ref_) + 3sol_std) .*100 
     push!(stds_matern3, sol_std)
     push!(errors_matern3, errors_)
end

reference_matern5 = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_matern5.t)
errors_matern5 = []
stds_matern5 = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_matern5.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_matern5.pu))[ch, :]
     sol_std = reduce(hcat, std.(sol_matern5.pu))[ch, :]
     errors_ = (abs.((sol_ .- ref_) ./ ref_) + 3sol_std) .*100 
     push!(stds_matern5, sol_std)
     push!(errors_matern5, errors_)
end



reference_ioup3 = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_ioup3.t)
errors_ioup3 = []
stds_ioup3 = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_ioup3.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_ioup3.pu))[ch, :]
     sol_std = reduce(hcat, std.(sol_ioup3.pu))[ch, :]
     errors_ = (abs.((sol_ .- ref_) ./ ref_) + 3sol_std) .*100 
     push!(stds_ioup3, sol_std)
     push!(errors_ioup3, errors_)
end

reference_ioup5 = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=sol_ioup5.t)
errors_ioup5 = []
stds_ioup5 = []
for ch in range(1,4) 
     ref_ = reduce(hcat, reference_ioup5.u)[ch, :]
     sol_ = reduce(hcat, mean.(sol_ioup5.pu))[ch, :]
     sol_std = reduce(hcat, std.(sol_ioup5.pu))[ch, :]
     errors_ = (abs.((sol_ .- ref_) ./ ref_) + 3sol_std) .*100 
     push!(stds_ioup5, sol_std)
     push!(errors_ioup5, errors_)
end


Plots.theme(
    :dao;
    markerstrokewidth=1.2,
    legend=:left,
)
colors = [1 1 2 2 3 3]
using Plots.PlotMeasures
p1 = plot(sol_matern3.t, errors_matern3[1], label="RosenbrockExpEK, \nMatern (3)", xlabel="t", ylabel="error (V), %", title="Absolute % error + 3std (V)", left_margin = [20mm 0mm], upper_margin = 40px, right_margin = [20mm 0mm])
p2 = plot(sol_matern3.t, errors_matern3[2], label="RosenbrockExpEK, \nMatern (3)", xlabel="t", ylabel="error (m), %", title="Error (m)", left_margin = [20mm 0mm], right_margin = [20mm 0mm])
p3 = plot(sol_matern3.t, errors_matern3[3], label="RosenbrockExpEK, \nMatern (3)", xlabel="t", ylabel="error (n), %", title="Error (n)", left_margin = [20mm 0mm], right_margin = [20mm 0mm])
p4 = plot(sol_matern3.t, errors_matern3[4], label="RosenbrockExpEK, \nMatern (3)", xlabel="t", ylabel="error (h), %", title="Error (h)", left_margin = [20mm 0mm], bottom_margin = 60px, right_margin = [20mm 0mm])

plot!(p1, sol_matern5.t, errors_matern5[1],label="EK1, Matern (3)", xlabel="t")
plot!(p1, sol_ioup3.t, errors_ioup3[1],label="RosenbrockExpEK, \nIOUP (3)", xlabel="t")
plot!(p1, sol_ioup5.t, errors_ioup5[1],label="EK1, IOUP (3)", xlabel="t", legend=false)
plot!(p1, sol_iwp3.t, errors_iwp3[1],label="RosenbrockExpEK, IWP (3)", xlabel="t", legend=false)
plot!(p1, sol_iwp5.t, errors_iwp5[1],label="EK1, IWP (3)", xlabel="t", legend=false)

plot!(p2, sol_matern5.t, errors_matern5[2],label="EK1, Matern (3)", xlabel="t")
plot!(p2, sol_ioup3.t, errors_ioup3[2],label="RosenbrockExpEK, \nIOUP (3)", xlabel="t")
plot!(p2, sol_ioup5.t, errors_ioup5[2],label="EK1, IOUP (3)", xlabel="t", legend=false)
plot!(p2, sol_iwp3.t, errors_iwp3[2],label="RosenbrockExpEK, IWP (3)", xlabel="t", legend=false)
plot!(p2, sol_iwp5.t, errors_iwp5[2],label="EK1, IWP (3)", xlabel="t", legend=false)

plot!(p3, sol_matern5.t, errors_matern5[3], xlabel="t")
plot!(p3, sol_ioup3.t, errors_ioup3[3], xlabel="t")
plot!(p3, sol_ioup5.t, errors_ioup5[3], xlabel="t", legend=false)
plot!(p3, sol_iwp3.t, errors_iwp3[3],label="RosenbrockExpEK, IWP (3)", xlabel="t", legend=false)
plot!(p3, sol_iwp5.t, errors_iwp5[3],label="EK1, IWP (3)", xlabel="t", legend=false)


plot!(p4, sol_matern5.t, errors_matern5[4],label="EK1, Matern (3)", xlabel="t")
plot!(p4, sol_ioup3.t, errors_ioup3[4],label="RosenbrockExpEK, \nIOUP (3)", xlabel="t")
plot!(p4, sol_ioup5.t, errors_ioup5[4],label="EK1, IOUP (4)", xlabel="t")
plot!(p4, sol_iwp3.t, errors_iwp3[4],label="RosenbrockExpEK, IWP (3)", xlabel="t")
plot!(p4, sol_iwp5.t, errors_iwp5[4],label="EK1, IWP (3)", xlabel="t")


l = @layout [a ; b ; c ; d]
plot!(p4, legendfont=font(7))
plot(p1, p2, p3, p4, layout = l)
plot!(size=(1000,600))
