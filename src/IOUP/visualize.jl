include("../utils.jl")
using Plots
using LinearAlgebra

evaltspan = (0.0, 5.0)
tstops = range(evaltspan..., length=500)

l1= -25
l2 = -0.1
l3= 0.1
l4 = 25 
n_derivatives = 3

names = [
        "IOUP($n_derivatives, $l1)"
        "IOUP($n_derivatives, $l2)"
        "IOUP($n_derivatives, $l3)"
        "IOUP($n_derivatives, $l4)"
        ]

priors = [
        IOUP(n_derivatives, l1),
        IOUP(n_derivatives, l2),
        IOUP(n_derivatives, l3),
        IOUP(n_derivatives, l4),
        ]


to_plot= []
for i in 1:length(priors)
    prior = priors[i]
    name = names[i]
    push!(to_plot, plot(prior, tstops; title=names[i], ylabel="", ylims=(-20,20)))
end

plot(to_plot...)
savefig("./visuals/IOUP/IOUP-priors-rates.png")



l1= 10^1.5
l2 = 10^0.85
l3 = 10^1.15
n_derivatives = 3



priors = [
        IOUP(n_derivatives, l1),
        IOUP(n_derivatives, l2),
        IOUP(n_derivatives, l3),
        IWP(n_derivatives),
        ]

l1 = round(Int64, l1)
l2 = round(Int64, l2)
l3 = round(Int64, l3)
names = [
        "IOUP($n_derivatives, $l1) - Equivalent loss to IWP"
        "IOUP($n_derivatives, $l2) -- min_1"
        "IOUP($n_derivatives, $l3) -- min_2"
        "IWP($n_derivatives)"
        ]

ylims_ = [
        (-10^10,10^10),
        (-10^8,10^8),
        (-10^10,10^10),
        (-10,10),
        ]
to_plot= []
for i in 1:length(priors)
    prior = priors[i]
    name = names[i]
    push!(to_plot, plot(prior, tstops; title=names[i], ylabel="", ylims=ylims_[i]))
end

plot(to_plot...)

savefig("./visuals/IOUP/IOUP-comparison-IWP.png")


# TODO: initialize with stationary distribution 
# TODO: make plots for 1st derivative

# derivative = 1
# dt= 0.01

# axes = Plots.plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))

# for i in 1:length(algorithms)

#     alg = algorithms[i] 
#     name = names[i]
#     color = colors[i]

#     prior = alg.prior
#     q = prior.num_derivatives

#     Ah, Qh = ProbNumDiffEq.discretize(prior, dt)
#     m_Qh = Matrix(Qh)
#     m_Ah = Matrix(Ah)
    

#     E0 = ProbNumDiffEq.projection(d, q)(0)
#     if q > 2
#         E1 = ProbNumDiffEq.projection(d, q)(1)
#         E2 = ProbNumDiffEq.projection(d, q)(2)
#         E = (E0, E1, E2)
#     elseif q > 1
#         E1 = ProbNumDiffEq.projection(d, q)(1)
#         E = (E0, E1)
#     else
#         E = Tuple(E0)
#     end


#     # Initial mean and covariance
#     initial = ProbNumDiffEq.initial_distribution(prior)

#     m0 = mean(initial)
#     C0 = Matrix(Diagonal(var(initial)))

#     tspan = (0.0, 50.0)
#     ts = collect(range(tspan..., length=5000))

#     Ys = sample_prior(m0, C0, ts, m_Ah, m_Qh, d, q)

#     titles = ["Samples from IOUP(l=-3)", "", "", ""]
#     for _ in 1:10
#         Ys = sample_prior(m0, C0, ts, m_Ah, m_Qh, d, q)
#         axes = plot_sample(ts, Ys, E, axes, titles, derivative, color, d; label="")
#     end
#     axes = plot_sample(ts, Ys, E, axes, titles, derivative, color,d; label=name)
    
# end

# display(axes)

# if TO_SAVE
#     Plots.savefig(axes, path_to_save)
# else
#     display(axes)
# end


