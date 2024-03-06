include("../../src/utils.jl")
using  LaTeXStrings

dts = 10.0 .^ range(-1, -4, length=11)[1:end-1]
evaltspan = (0.0, 50.0)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()

_setups = [
    "Tsit5" => Dict(:alg => Tsit5())
    "RK4" => Dict(:alg=>RK4())
    "Implicit Euler" => Dict(:alg=>ImplicitEuler())
    "Exponential Euler" => Dict(:alg=>LawsonEuler(krylov = true, m = 50)) # on krylov and m : https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/ 
    "EK1(3)" => Dict(:alg=>EK1(order=3, smooth=DENSE, diffusionmodel=DM))
    ]
setups = last.(_setups)

# COLORS
colors = []

# c1 = colorant"orange"
# c2 = colorant"green"
# colors2 = range(c1, stop=c2, length=4)

# for c in colors2
#     push!(colors, c)
# end
push!(colors, colorant"orange")
push!(colors, colorant"tomato1")
push!(colors, colorant"green")
push!(colors, colorant"deepskyblue")
push!(colors, colorant"black")

colors

labels = [
    "Tsit5",
    "RK4", 
    "Implicit Euler",
    "Exponential Euler",
    "Probabilistic solver",
    "", 
    "",
    "",
    "",
]


prob = define_problem()
algorithms = [setup[:alg] for setup in setups]

p = plot(legendfont = font(7), size = (500, 400))

for i in 1:length(algorithms)
    algorithm = algorithms[i]
    errors = []

    for dt in dts
        try
            solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)
            reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
            error = l2(solution.u, reference(solution.t))
            push!(errors, error)
        catch e
            push!(errors, NaN)
            continue
        end
    
    end
    
    plot!(p, dts, errors, label=labels[i], 
        framestyle=:axes,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)

    scatter!(p, dts, errors, label="", 
        framestyle=:axes,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)
    display(p)
end



prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()

for i in 1:5
    algorithm = algorithms[i]
    errors = []

    for dt in dts
        try
            solution_transformed = solve(prob_transformed,  algorithm, dt=dt, adaptive=false, dense = false)
            reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
            error = l2_transformed(solution_transformed, reference, inverse_transforms)
            push!(errors, error)    
        catch e
            push!(errors, NaN)
            continue
        end
    
    end
    if i == 5
        transformed_label = "Transformed"
    else
        transformed_label = ""
    end
    plot!(p, dts, errors, label=transformed_label, 
        framestyle=:axes,
        linestyle=:dash,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)

    scatter!(p, dts, errors, label="", 
        framestyle=:axes,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)
    display(p)
end

plot!(p, xaxis=:log10, yaxis=:log10, 
        legend=:bottomright, xlabel=L"$\Delta$t", 
        ylabel="tRMSE", dpi=600)
p

savefig(p, "./visuals/deterministic/deterministic-transformation.pdf")