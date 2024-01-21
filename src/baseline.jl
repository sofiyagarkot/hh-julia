# experiments/January/week-2/fixed_diffusion.jl:
using DifferentialEquations
# using GeometricIntegratorsDiffEq
include("utils.jl")

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = true

# fixed diffusion
DM = FixedDiffusion()

# EK1  comparing to Exponential Euler, and (optionally) BackwardEuler

c1 = colorant"lightblue"
c2 = colorant"darkblue"
colors = range(c1, stop=c2, length=4)
push!(colors, colorant"darkorange")

algorithms = [
             EK1(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
             EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
             LawsonEuler(krylov = true, m = 50),   # https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/#Semilinear-ODE 
            #  GIImplicitEuler(),
             ]


names = [  
        "IWP(1)", 
        "IWP(2)", 
        "IWP(3)", 
        "IWP(4)",
        "Exponential Euler",
        # "Implicit Euler",
        ]


prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

solutions = get_solutions(algorithms[1:length(algorithms)-1], prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE) 
solution_euler = solve(prob, algorithms[end], dt=0.01, saveat=tstops)
push!(solutions, solution_euler)

errors, stds = absolute_errors(solutions, prob)
errors_prob = errors[1:length(errors)-1]
error_euler =  errors[length(errors)]
times = [sol_i.t for sol_i in solutions[1:length(errors)-1]]

# Absolute errors
titles=["Absolute error with confidence intervals (EK1, FixedDiffusion)", "", "", ""];

p = plot_errors(
    [solution_euler.t], [error_euler], 
    ["Exponential Euler"], titles;
    colors=[colors[end]],
    to_save_path="./visuals/baseline/absolute_errors_in_time.png",
    to_save=TO_SAVE)
   
plot_errors(
    times, errors_prob, names[1:length(algorithms)-1], titles; 
    stds, 
    to_save_path="./visuals/baseline/absolute_errors_in_time.png", 
    to_save=TO_SAVE,
    p = p, 
    colors=colors[1:length(errors)-1]
    )


# Log absolute errors
titles=["Log of absolute errors (EK1, dt=0.01, FixedDiffusion)", "", "", ""];

p = plot_errors(
        times, errors_prob, names[1:length(algorithms)-1], titles; 
        stds, 
        to_save_path="./visuals/baseline/log_absolute_errors_in_time.png", 
        to_save=TO_SAVE,
        log_plot = true,
        colors=colors[1:length(errors)-1]
        )

plot_errors(
    [solution_euler.t], [error_euler], 
    ["Exponential Euler"], titles;
    colors=[colors[end]],
    to_save_path="./visuals/baseline/log_absolute_errors_in_time.png",
    p = p, 
    log_plot = true,
    to_save=TO_SAVE)


# Work-Precision plots

dts = 10.0 .^ range(-2, -7, length=11)[begin:end-1]
abstols = reltols = repeat([missing], length(dts))


plots = work_precision_plot(
        names, algorithms, prob; 
        DENSE = DENSE, 
        SAVE_EVERYSTEP = false, 
        to_save=TO_SAVE, 
        to_save_path="./visuals/baseline/fixed_diffusion_wp_EK1_IWP.png",
        to_save_path2 = "./visuals/baseline/fixed_diffusion_steps_number_wp_EK1_IWP.png",
        title="EK1, dt = 0.01, FixedDiffusion",
        adaptive=ADAPTIVE,
        abstols=abstols,
        reltols=reltols,
        dts=dts,
        colors = colors
        )

using Plots
to_save_path="./visuals/baseline/fixed_diffusion_wp_EK1_IWP.png"
to_save_path2 = "./visuals/baseline/fixed_diffusion_steps_number_wp_EK1_IWP.png"
Plots.savefig(plots[1], to_save_path)
Plots.savefig(plots[2], to_save_path2)
