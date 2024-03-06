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

numruns = 20
time_tmp = Vector{Float64}(undef, numruns)

algorithms = [setup[:alg] for setup in setups]

p = plot(legendfont = font(7), size = (500, 400))

for i in 1:length(algorithms)
    algorithm = algorithms[i]
    errors = []
    times = []

    for dt in dts
        try
            solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)
            reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
            error = l2(solution.u, reference(solution.t))
            push!(errors, error)

            for i in 1:numruns
                time_tmp[i] =  @elapsed solve(prob_transformed, algorithm, dt=dt, adaptive=false, dense = false)
            end
            push!(times, mean(time_tmp))

        catch e
            push!(errors, NaN)
            push!(times, NaN)
            continue           
        end
    
    end
    
    plot!(p, errors, times, label=labels[i], 
        framestyle=:axes,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)

    scatter!(p, errors, times, label="", 
        framestyle=:axes,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)
    display(p)
end

p
prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()

for i in 1:length(algorithms)
    algorithm = algorithms[i]
    errors = []
    times = []

    for dt in dts
        try
            solution_transformed = solve(prob_transformed,  algorithm, dt=dt, adaptive=false, dense = false)
            reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_transformed.t)
            error = l2_transformed(solution_transformed, reference, inverse_transforms)
            push!(errors, error) 
            
            
            for i in 1:numruns
                time_tmp[i] =  @elapsed solve(prob_transformed, algorithm, dt=dt, adaptive=false, dense = false)
            end
            push!(times, mean(time_tmp))

        catch e
            push!(errors, NaN)
            push!(times, NaN)
            continue
        end
    
    end
    if i == 2
        transformed_label = "Transformed"
    else
        transformed_label = ""
    end
    plot!(p, errors, times, label=transformed_label, 
        framestyle=:axes,
        linestyle=:dash,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)

    scatter!(p, errors, times, label="", 
        framestyle=:axes,
        xaxis=:log10, yaxis=:log10,
        color=colors[i],
        fg_legend = :transparent)

    display(p)
end

plot!(p, xaxis=:log10, yaxis=:log10, 
        legend=:topright, legendtitle="Prior of probabilistic solver",
        legendtitlefontsize=8,
        legendfontsize=7,
        ylabel="Time, s ", 
        xlabel="tRMSE", dpi=600)
p

savefig(p, "./visuals/iwp/transformed-time-error.pdf")

