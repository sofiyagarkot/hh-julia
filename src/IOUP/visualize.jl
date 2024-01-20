include("../utils.jl")
using Plots
using LinearAlgebra

# TODO: initialize with static distribution, 
# TODO: make plots for several rate params as 4x4 grid


SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = true


prob = define_problem()

l = -3
d = 4     # dimension of ODE , V, m , n , h   
names = [
        "IOUP(1, $l)"
        "IOUP(2, $l)"
        "IOUP(3, $l)"
        "IOUP(4, $l)"
        ]

algorithms = [
    # EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false) 
    EK1(prior=IOUP(d, 1, l,)),
    EK1(prior=IOUP(d, 2, l,)),
    EK1(prior=IOUP(d, 3, l,)),
    EK1(prior=IOUP(d, 4, l,)),
    ]

dt = 0.01
derivative = 1
uElType = Float64
path_to_save = "./visuals/January/plots_of_priors/IOUP/$derivative-derivative-l-$l-orders.png"
colors = (1:length(algorithms))


axes = Plots.plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))
Qhs = [] 
for i in 1:length(algorithms)

    alg = algorithms[i] 
    name = names[i]
    color = colors[i]

    prior = alg.prior
    q = prior.num_derivatives

    Ah, Qh = ProbNumDiffEq.discretize(prior, dt)
    m_Qh = Matrix(Qh)
    m_Ah = Matrix(Ah)
    

    E0 = ProbNumDiffEq.projection(d, q)(0)
    if q > 2
        E1 = ProbNumDiffEq.projection(d, q)(1)
        E2 = ProbNumDiffEq.projection(d, q)(2)
        E = (E0, E1, E2)
    elseif q > 1
        E1 = ProbNumDiffEq.projection(d, q)(1)
        E = (E0, E1)
    else
        E = Tuple(E0)
    end


    # Initial mean and covariance
    m0 = zeros(d * (q + 1))
    C0 = float(Matrix(I, d * (q + 1), d * (q + 1)))

    tspan = (0.0, 5.0)
    ts = collect(range(tspan..., length=500))

    Ys = sample_prior(m0, C0, ts, m_Ah, m_Qh, d, q)

    titles = ["Samples from IOUP(l=-3)", "", "", ""]
    for _ in 1:10
        Ys = sample_prior(m0, C0, ts, m_Ah, m_Qh, d, q)
        axes = plot_sample(ts, Ys, E, axes, titles, derivative, color; label="")
    end
    axes = plot_sample(ts, Ys, E, axes, titles, derivative, color; label=name)
    
end

display(axes)
TO_SAVE = true
if TO_SAVE
    Plots.savefig(axes, path_to_save)
else
    display(axes)
end


