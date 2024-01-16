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

function get_drift_and_diffusion_matrices(p::ProbNumDiffEq.Matern)
    q = p.num_derivatives
    l = p.lengthscale
    ν = q - 1 / 2
    λ = sqrt(2ν) / l

    F_breve = diagm(1 => ones(q))
    @. F_breve[end, :] = -binomial(q + 1, 0:q) * λ^((q+1):-1:1)

    L_breve = zeros(q + 1)
    L_breve[end] = 1.0

    d = p.wiener_process_dimension
    # F = ProbNumDiffEq.IsometricKroneckerProduct(d, F_breve)
    # L = ProbNumDiffEq.IsometricKroneckerProduct(d, L_breve)
    F = kron(I(d), F_breve)
    L = kron(I(d), L_breve)

    return F, L
end

function get_drift_and_diffusion_matrices(p::ProbNumDiffEq.IWP)
    d, q = p.wiener_process_dimension, p.num_derivatives

    F_breve = diagm(1 => ones(q))
    L_breve = zeros(q + 1)
    L_breve[end] = 1.0

    # F = ProbNumDiffEq.IsometricKroneckerProduct(d, F_breve)
    # L = ProbNumDiffEq.IsometricKroneckerProduct(d, L_breve)
    F = kron(I(d), F_breve)
    L = kron(I(d), L_breve)
    return F, L
end

function get_drift_and_diffusion_matrices(p::ProbNumDiffEq.IOUP)
    q = p.num_derivatives
    r = p.rate_parameter

    F_breve = diagm(1 => ones(q))
    F_breve[end, end] = r

    L_breve = zeros(q + 1)
    L_breve[end] = 1.0

    d = p.wiener_process_dimension
    # F = ProbNumDiffEq.IsometricKroneckerProduct(d, F_breve)
    # L = ProbNumDiffEq.IsometricKroneckerProduct(d, L_breve)
    F = kron(I(d), F_breve)
    L = kron(I(d), L_breve)
    return F, L
end