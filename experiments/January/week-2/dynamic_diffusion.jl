include("../../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = false

# dynamic diffusion

# EK1 
algorithms = [
             EK1(prior=IWP(1), smooth=SMOOTH),
             EK1(prior=IWP(2), smooth=SMOOTH),
             EK1(prior=IWP(3), smooth=SMOOTH), 
             EK1(prior=IWP(4), smooth=SMOOTH), 
             ]

names = [  
        "IWP(1)", 
        "IWP(2)", 
        "IWP(3)", 
        "IWP(4)", 
        # "IWP(5)"
        ]


prob = define_problem()
solutions = get_solutions(algorithms, prob; adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

TO_SAVE = true
titles=["Absolute errors with confidence intervals (EK1, adaptive step sizes, DynamicDiffusion)", "", "", ""];

errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)

plot_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/absolute_errors/dynamic_diffusions_errors_EK1.png", 
        to_save=TO_SAVE
        )
        
titles=["Log of absolute errors with confidence intervals (EK1, adaptive step sizes, DynamicDiffusion)", "", "", ""];

plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/absolute_errors/log_dynamic_diffusions_errors_EK1_ribbon.png", 
        to_save=TO_SAVE
        )


plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/absolute_errors/log_dynamic_diffusions_absolute_errors_EK1.png", 
        to_save=TO_SAVE,
        ribbon = false
        )

# work_precision_plot(
#         names, algorithms, prob; 
#         DENSE = DENSE, 
#         SAVE_EVERYSTEP = false, 
#         to_save=TO_SAVE, 
#         to_save_path="./visuals/January/diffusions/dynamic_diffusion_wp_EK1_IWP.png",
#         title="EK1, DynamicDiffusion",
#         adaptive=ADAPTIVE,
#         )


# EK0 
algorithms = [
    EK0(prior=IWP(1), smooth=SMOOTH),
    EK0(prior=IWP(2), smooth=SMOOTH),
    EK0(prior=IWP(3), smooth=SMOOTH), 
    EK0(prior=IWP(4), smooth=SMOOTH), 
#      EK1(prior=IWP(5), smooth=SMOOTH, diffusionmodel=DM)
    ]

names = [  
"IWP(1)", 
"IWP(2)", 
"IWP(3)", 
"IWP(4)", 
# "IWP(5)"
]


prob = define_problem()

solutions = get_solutions(algorithms, prob; adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

TO_SAVE = true
titles=["Absolute errors with confidence intervals (EK0, adaptive, DynamicDiffusion)", "", "", ""];

plot_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/dynamic_diffusions_errors_EK0.png", 
    to_save=TO_SAVE
    )

# work_precision_plot(
#     names, algorithms, prob; 
#     DENSE = DENSE, 
#     SAVE_EVERYSTEP = false, 
#     to_save=TO_SAVE, 
#     to_save_path="./visuals/January/diffusions/absolute_errors/dynamic_diffusion_wp_EK0_IWP.png",
#     title="EK0, DynamicDiffusion, adaptive",
#     adaptive=ADAPTIVE,
#     tstops=tstops
#     )

titles=["Log of absolute % error with confidence intervals (EK0, adaptive, DynamicDiffusion)", "", "", ""];

plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/absolute_errors/log_dynamic_diffusions_errors_EK0_ribbon.png", 
        to_save=TO_SAVE
        )


plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/log_dynamic_diffusions_absolute_errors_EK0.png", 
    to_save=TO_SAVE,
    ribbon = false
    )

# Comparing EK0 and EK1 errors

algorithms = [
    EK0(prior=IWP(2), smooth=SMOOTH),
    EK0(prior=IWP(3), smooth=SMOOTH),
    EK0(prior=IWP(4), smooth=SMOOTH),
    EK0(prior=IWP(5), smooth=SMOOTH),
    EK1(prior=IWP(2), smooth=SMOOTH), 
    EK1(prior=IWP(3), smooth=SMOOTH), 
    EK1(prior=IWP(4), smooth=SMOOTH), 
    EK1(prior=IWP(5), smooth=SMOOTH), 
    ]

names = [  
    "EK0(2)", 
    "EK0(3)",  
    "EK0(4)",     
    "EK0(5)", 
    "EK1(2)", 
    "EK1(3)", 
    "EK1(4)", 
    "EK1(5)", 
    ]

work_precision_plot(
    names, algorithms, prob; 
    DENSE = DENSE, 
    SAVE_EVERYSTEP = false, 
    to_save=TO_SAVE, 
    to_save_path="./visuals/January/diffusions/dynamic_diffusion_wp_EK0_vs_EK1.png",
    title="EK1 vs. EK0, DynamicDiffusion, adaptive",
    adaptive=ADAPTIVE,
    tstops=tstops,
    colors = [:red :red :red :red :blue :blue :blue :blue],
    )


algorithms = [
    # EK0(prior=IWP(2), smooth=SMOOTH),
    EK0(prior=IWP(3), smooth=SMOOTH),
    # EK0(prior=IWP(4), smooth=SMOOTH),
    # EK0(prior=IWP(5), smooth=SMOOTH),
    # EK1(prior=IWP(2), smooth=SMOOTH), 
    EK1(prior=IWP(3), smooth=SMOOTH), 
    # EK1(prior=IWP(4), smooth=SMOOTH), 
    # EK1(prior=IWP(5), smooth=SMOOTH), 
    ]

names = [  
    # "EK0(2)", 
    "EK0(3)",  
    # "EK0(4)",     
    # "EK0(5)", 
    # "EK1(2)", 
    "EK1(3)", 
    # "EK1(4)", 
    # "EK1(5)", 
    ]
titles = ["Log abs error (DynamicDiffusion, adaptive step sizes)", "", "", ""]


solutions = get_solutions(algorithms, prob; adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/diffusions/absolute_errors/log_dynamic_diffusions_errors_EK0-EK1_ribbon.png", 
    to_save=false
    )



ch = 1
evaltspan = (0.0, 50.1)
ts = range(evaltspan..., length=5000)
derivative = 1
# plot_solution_with_uncertainty( solutions[2], ch, 0, "EK0(3) Voltage", ts=ts, to_save =false, to_save_path="./visuals/January/diffusions/solution-fixed_diffusion_EK0(3)_voltage.png")
# plot_solution_with_uncertainty( solutions[5], ch, 0, "EK1(3) Voltage", ts=ts, to_save =false, to_save_path="./visuals/January/diffusions/solution-fixed_diffusion_EK1(3)_voltage.png")
titles= ["Comparing the solutions for dynamic diffusion,"*" $derivative derivative", "", "", ""]
plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        to_save=true,
        to_save_path="./visuals/January/diffusions/solution-dynamic_diffusion_1st_derivative_EK0-EK1_voltage.png"
        )
        
ch = 1
evaltspan = (0.0, 56.1)
ts = range(evaltspan..., length=5000)
derivative = 1
titles = ["Log std for dynamic diffusion,"*" $derivative derivative", "", "", ""]
plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        log_std=true,
        to_save=true,
        to_save_path="./visuals/January/diffusions/log_std-dynamic_diffusion_1st_derivative_EK0-EK1_voltage.png"
        )
