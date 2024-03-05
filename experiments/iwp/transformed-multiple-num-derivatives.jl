using RecursiveArrayTools
include("../../src/utils.jl")
using  LaTeXStrings

dts = 10.0 .^ range(-1, -4, length=11)[1:end-1]
evaltspan = (0.0, 50.0)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()

_setups = [
    "EK1(2)" => Dict(:alg=>EK1(order=2, smooth=DENSE, diffusionmodel=DM))
    "EK1(3)" => Dict(:alg=>EK1(order=3, smooth=DENSE, diffusionmodel=DM))
    "EK1(4)" => Dict(:alg=>EK1(order=4, smooth=DENSE, diffusionmodel=DM))
    ]
    
setups = last.(_setups)

# COLORS
colors = ["blue2", "red3", "coral"]


labels = [
    "IWP(2)",
    "IWP(3)",
    "IWP(4)",
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
        catch e
            push!(errors, NaN)
            continue

        else            
            solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)
            reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
            error = l2(solution.u, reference(solution.t))
            push!(errors, error)
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
end



prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()

for i in 1:length(algorithms)
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
    if i == 2
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
end

plot!(p, xaxis=:log10, yaxis=:log10, 
        legend=:bottomright, legendtitle="Prior of probabilistic solver",
        legendtitlefontsize=8,
        legendfontsize=7,
         xlabel=L"$\Delta$t", 
        ylabel="tRMSE", dpi=600)
p

savefig(p, "./visuals/iwp/transformed-multiple-num-derivatives.png")