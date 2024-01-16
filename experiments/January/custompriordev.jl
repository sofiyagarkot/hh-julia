using ProbNumDiffEq, LinearAlgebra, Plots

"""
Linear time-invariant SDE:
dX(t) = F X(t) dt + L dW(t)
"""
struct CustomSDEPrior{elType,FT,LT} <: ProbNumDiffEq.AbstractODEFilterPrior{elType}
    # These two fields are required for the prior to work
    wiener_process_dimension::Int
    num_derivatives::Int
    # This is the custom stuff
    F::FT
    L::LT
end
# Constructors
CustomSDEPrior(wiener_process_dimension, num_derivatives, F, L) =
    CustomSDEPrior{Float64,typeof(F),typeof(L)}(wiener_process_dimension, num_derivatives, F, L)
# Functions that need to be implemented
ProbNumDiffEq.to_sde(p::CustomSDEPrior) = ProbNumDiffEq.LTISDE(p.F, p.L)
ProbNumDiffEq.discretize(p::CustomSDEPrior, dt::Real) =
    ProbNumDiffEq.discretize(ProbNumDiffEq.to_sde(p), dt)

# Try it out
function lotka_volterra(du, u, p, t)
    du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = -p[3] * u[2] + p[4] * u[1] * u[2]
end
p = [1.5, 1.0, 3.0, 1.0]
u0 = [1.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)

sol1 = solve(prob, EK1(order=1));
# plot(sol1)

# Specify an integrated wiener process myself:
F, L = kron(I(2), [0.0 1.0; 0.0 0.0]), kron(I(2), [0.0; 1.0])
sol2 = solve(prob, EK1(prior=CustomSDEPrior(2, 1, F, L)));

sol1[end] ≈ sol2[end]
