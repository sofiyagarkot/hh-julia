include("../../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = true

# fixed diffusion
DM = FixedDiffusion()

# EK1 
algorithms = [
             EK1(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
             EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
             ]

names = [  
        "IWP(1)", 
        "IWP(2)", 
        "IWP(3)", 
        "IWP(4)", 
        ]


prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

titles=["Absolute error with confidence intervals (EK1, FixedDiffusion)", "", "", ""];

plot_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/absolute_errors/fixed_diffusions_errors_EK1.png", 
        to_save=TO_SAVE
        )
        

dts = 10.0 .^ range(-2, -5, length=10)[begin:end-1]
abstols = reltols = repeat([missing], length(dts))
work_precision_plot(
        names, algorithms, prob; 
        DENSE = DENSE, 
        SAVE_EVERYSTEP = false, 
        to_save=TO_SAVE, 
        to_save_path="./visuals/January/diffusions/fixed_diffusion_wp_EK1_IWP.png",
        title="EK1, FixedDiffusion",
        adaptive=ADAPTIVE,
        abstols=abstols,
        reltols=reltols,
        dts=dts,
        colors = [1 1 1 2 2 2]
        )
        
titles=["Log of absolute errors with confidence intervals (EK1, dt=0.01, FixedDiffusion)", "", "", ""];

plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/absolute_errors/log_fixed_diffusions_errors_EK1_ribbon.png", 
        to_save=TO_SAVE
        )

plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/log_fixed_diffusions_errors_EK1.png", 
    to_save=TO_SAVE,
    ribbon = false
    )

# EK0 
algorithms = [
    EK0(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
    EK0(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
    ]

names = [  
"IWP(1)", 
"IWP(2)", 
"IWP(3)", 
"IWP(4)", 
]


prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

titles=["Absolute error with confidence intervals (EK0, dt=0.01, FixedDiffusion)", "", "", ""];

plot_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/fixed_diffusions_errors_EK0.png", 
    to_save=TO_SAVE
    )

work_precision_plot(
        names, algorithms, prob; 
        DENSE = DENSE, 
        SAVE_EVERYSTEP = false, 
        to_save=TO_SAVE, 
        to_save_path="./visuals/January/diffusions/fixed_diffusion_wp_EK0_IWP.png",
        title="EK0, FixedDiffusion",
        adaptive=ADAPTIVE,
        abstols=abstols,
        reltols=reltols,
        dts=dts
        )

titles=["Log of absolute error with confidence intervals (EK0, dt=0.01, FixedDiffusion)", "", "", ""];

plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/log_fixed_diffusions_errors_EK0_ribbon.png", 
    to_save=TO_SAVE
    )


plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/log_fixed_diffusions_errors_EK0.png", 
    to_save=TO_SAVE,
    ribbon = false
    )

# Comparing EK0 and EK1 errors
algorithms = [
    EK0(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM),
    EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM), 
    EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
    EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
    ]

names = [  
    "EK0(2)", 
    "EK0(3)",  
    "EK0(4)",  
    "EK1(2)", 
    "EK1(3)", 
    "EK1(4)", 
    ]


work_precision_plot(
    names, algorithms, prob; 
    DENSE = DENSE, 
    SAVE_EVERYSTEP = false, 
    to_save=TO_SAVE, 
    to_save_path="./visuals/January/diffusions/fixed_diffusion_wp_EK0_vs_EK1_IWP.png",
    title="EK0 vs EK1, FixedDiffusion",
    adaptive=ADAPTIVE,
    abstols=abstols,
    reltols=reltols,
    dts=dts
    )

# dts = exp10.(-1:-1:-5)
# plot_errors_vs_dt(dts, prob, algorithms, names ; 
#                 TO_SAVE=TO_SAVE, 
#                 to_save_path="./visuals/January/diffusions/fixed_diffusion_errors_vs_dt.png",
#                 dense=DENSE,
#                 )
#  _________________________________________________________________

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = true

# fixed diffusion
DM = FixedDiffusion()

# Comparing EK0 and EK1 errors
algorithms = [
    # EK0(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
    EK0(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM),
    # EK0(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM),
    # EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM), 
    EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
    # EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
    ]

names = [  
    # "EK0(2)", 
    "EK0(3)",  
    # "EK0(4)",  
    # "EK1(2)", 
    "EK1(3)", 
    # "EK1(4)", 
    ]

evaltspan = (0.0, 50.0)
dts = 10.0 .^ range(-2, -3, length=10)[begin:end-1]
abstols = reltols = repeat([missing], length(dts))

evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

# work_precision_plot(
#     names, algorithms, prob; 
#     DENSE = DENSE, 
#     SAVE_EVERYSTEP = false, 
#     to_save=TO_SAVE, 
#     to_save_path="./visuals/January/diffusions/fixed_diffusion_wp_EK0_vs_EK1_IWP.png",
#     title="EK0 vs EK1, FixedDiffusion",
#     adaptive=ADAPTIVE,
#     abstols=abstols,
#     reltols=reltols,
#     dts=dts,
#     colors = [1 1 1 2 2 2]
#     )

titles = ["Log absolute error with log uncertainty (1 std)", "", "", ""]
plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/log_fixed_diffusions_errors_EK0-EK1_ribbon.png", 
    to_save=true
    )


ch = 1
evaltspan = (0.0, 50.1)
ts = range(evaltspan..., length=5000)
derivative = 1
# plot_solution_with_uncertainty( solutions[2], ch, 0, "EK0(3) Voltage", ts=ts, to_save =false, to_save_path="./visuals/January/diffusions/solution-fixed_diffusion_EK0(3)_voltage.png")
# plot_solution_with_uncertainty( solutions[5], ch, 0, "EK1(3) Voltage", ts=ts, to_save =false, to_save_path="./visuals/January/diffusions/solution-fixed_diffusion_EK1(3)_voltage.png")
titles= ["Comparing the solutions,"*" $derivative derivative", "", "", ""]
plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
     ts = ts,
     to_save=true,
     to_save_path="./visuals/January/diffusions/solution-fixed_diffusion_1st_derivative_EK0-EK1_voltage.png"
     )
     
ch = 1
evaltspan = (0.0, 56.1)
ts = range(evaltspan..., length=5000)
derivative = 1
titles = ["Log std for Fixed Diffusion,"*" $derivative derivative", "", "", ""]
plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        log_std=true,
        to_save=true,
        to_save_path="./visuals/January/diffusions/log_std-fixed_diffusion_1st_EK0-EK1_voltage.png"
        )
