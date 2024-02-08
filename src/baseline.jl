using RecursiveArrayTools
include("utils.jl")
using  LaTeXStrings

dts = 10.0 .^ range(-1, -5, length=11)[1:end-1]
evaltspan = (0.0, 50.0)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()

_setups = [
    "Tsit5" => Dict(:alg => Tsit5())
    # "Vern7" => Dict(:alg => Vern7())
    "RK4" => Dict(:alg=>RK4())
    "Implicit Euler" => Dict(:alg=>ImplicitEuler())
    "Exponential Euler" => Dict(:alg=>LawsonEuler(krylov = true, m = 50)) # on krylov and m : https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/ 
    "EK1(1)" => Dict(:alg=>EK1(order=1, smooth=DENSE, diffusionmodel=DM))
    "EK1(2)" => Dict(:alg=>EK1(order=2, smooth=DENSE, diffusionmodel=DM))
    "EK1(3)" => Dict(:alg=>EK1(order=3, smooth=DENSE, diffusionmodel=DM))
    "EK1(4)" => Dict(:alg=>EK1(order=4, smooth=DENSE, diffusionmodel=DM))
    ]
    
labels = first.(_setups)
setups = last.(_setups)

# COLORS
colors = [colorant"yellow"]
push!(colors, colorant"darkorange")
push!(colors, colorant"red2")
push!(colors, colorant"darkred")

c1 = colorant"lightblue"
c2 = colorant"darkblue"
colors2 = range(c1, stop=c2, length=4)

for c in colors2
    push!(colors, c)
end

# L2 loss from https://github.com/SciML/DiffEqDevTools.jl/blob/master/src/test_solution.jl 
function l2(sol, analytic)
    return sqrt(recursive_mean(vecvecapply((x) -> float(x) .^ 2,
                                                      sol - analytic)))
end

prob = define_problem()


algorithms = [setup[:alg] for setup in setups]

p = plot(legendfont = font(7), size = (500, 400))

for i in 1:length(algorithms)
    print("i: ", i, "\n")
    algorithm = algorithms[i]
    errors = []

    for dt in dts
        print("dt", dt, "\n")
        # tstops = range(evaltspan..., length=Int(round(evaltspan[end]/dt, digits = 0)))
        # print("tstops: ", tstops, "\n")
        try
            solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)
        catch e
            print("error: ", e, "\n")
            push!(errors, NaN)
            continue

        else            
            solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)

            # print("solution: ", solution, "\n")

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
    display(p)
end

plot!(p, xaxis=:log10, yaxis=:log10, 
        legend=:bottomright, xlabel=L"$\Delta$t", 
        ylabel="L2 error")
p

# savefig(p, "./visuals/my-wp-diagram.png")