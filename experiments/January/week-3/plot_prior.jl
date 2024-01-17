include("../../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = false


prob = define_problem()

names = [
        # "IOUP(3, 1), update=false",
        "IWP(3)"
        ]

algorithms = [
    # EK1(prior=IOUP(3,  update_rate_parameter=true)),
    EK1(prior=IWP(4,3), smooth=SMOOTH),
    ]

derivative = 0
dt = 0.01

uElType = Float64
# d = 4     # dimension of ODE , V, m , n , h

for alg in algorithms 
    q = 3
    # d = is_secondorder_ode ? length(u[1, :]) : length(u)
    d = 4     # dimension of ODE , V, m , n , h    
    FAC = ProbNumDiffEq.get_covariance_structure(alg; elType=Float64, d, q)

    # A, Q, Ah, Qh, P, PI = ProbNumDiffEq.initialize_transition_matrices(FAC, prior, dt)
    A, Q, Ah, Qh, P, PI = ProbNumDiffEq.initialize_transition_matrices(
        ProbNumDiffEq.DenseCovariance{Float64}(d, q), prior, dt)

    E0 = kron([1 0 0], Matrix(I, d, d))
    E1 = kron([0 1 0], Matrix(I, d, d))
    E2 = kron([0 0 1], Matrix(I, d, d))


    # Initial mean and covariance
    m0 = zeros(d * (q + 1))
    C0 = Matrix(I, d * (q + 1), d * (q + 1))

    tspan = (0.0, 5.0)
    ts = range(tspan..., length=500)
    
    fig, axes = plot_sample(ts, Ys)
    for _ in 1:10
        Ys = sample_prior(m0, C0, ts, A, Q)
        plot_sample(ts, Ys, axes=axes)
    end

end

plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        log_std=false,
        to_save=TO_SAVE,
        to_save_path="./visuals/January/ioup/solutions-fixed_diffusion_1st_voltage.png"
        )

ch = 1
evaltspan = (0.0, 5.0)
ts = range(evaltspan..., length=500)
derivative = 0
titles = ["Solutions for"*" $derivative derivative", "", "", ""]

plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        log_std=false,
        to_save=TO_SAVE,
        to_save_path="./visuals/January/ioup/solutions-fixed_diffusion_1st_voltage.png"
        )