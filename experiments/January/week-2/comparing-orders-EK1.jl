include("../../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = true

# fixed diffusion
DM = FixedDiffusion()

algorithms = [
             EK1(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
             EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
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
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=100000)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

TO_SAVE = true
titles=["Absolute % error with confidence intervals (EK1, dt=0.005, FixedDiffusion)", "", "", ""];

plot_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/fixed_diffusions_errors_EK1.png", 
        to_save=TO_SAVE
        )

work_precision_plot(
        names, algorithms, prob; 
        DENSE = DENSE, 
        SAVE_EVERYSTEP = false, 
        to_save=TO_SAVE, 
        to_save_path="./visuals/January/diffusions/fixed_diffusion_wp_EK1_IWP.png",
        title="EK1, FixedDiffusion",
        adaptive=ADAPTIVE,
        tstops=tstops
        )
        
titles=["Log of absolute % error with confidence intervals (EK1, dt=0.005, FixedDiffusion)", "", "", ""];

plot_log_errors(
        times, errors, names, titles; 
        # stds, 
        to_save_path="./visuals/January/diffusions/log_fixed_diffusions_errors_EK1.png", 
        to_save=TO_SAVE
        )

# dynamic diffusion

titles=["Absolute % error with confidence intervals (EK1, dt=0.005, DynamicDiffusion)", "", "", ""];
names = [  
        "IWP(1)", 
        "IWP(2)", 
        "IWP(3)", 
        "IWP(4)", 
        # "IWP(5)"
        ]
algorithms = [
                EK1(prior=IWP(1), smooth=SMOOTH),
                EK1(prior=IWP(2), smooth=SMOOTH),
                EK1(prior=IWP(3), smooth=SMOOTH), 
                EK1(prior=IWP(4), smooth=SMOOTH), 
                # EK1(prior=IWP(5), smooth=SMOOTH)
                ]


prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

plot_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/diffusions/dynamic_diffusions_errors.png", 
        to_save=TO_SAVE
        )

work_precision_plot(
        names, algorithms, prob; 
        DENSE = DENSE, 
        SAVE_EVERYSTEP = false, to_save=TO_SAVE, 
        to_save_path="./visuals/January/diffusions/dynamic_diffusion_wp_diagram-comparing-the-orders.png",
        title="EK1, dynamic diffusion"
        )


# on the same plot
names = [  
        "IWP(1)", 
        "IWP(1), fixed", 
        "IWP(2)",  
        "IWP(2), fixed", 
        "IWP(3)", 
        "IWP(3), fixed", 
        "IWP(4)", 
        "IWP(4), fixed", 
        ]
algorithms = [
             EK1(prior=IWP(1), smooth=SMOOTH),
             EK1(prior=IWP(1), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(2), smooth=SMOOTH),
             EK1(prior=IWP(2), smooth=SMOOTH, diffusionmodel=DM),
             EK1(prior=IWP(3), smooth=SMOOTH), 
             EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM), 
             EK1(prior=IWP(4), smooth=SMOOTH), 
             EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=DM), 
             ]


prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)  

std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
TO_SAVE = true
titles=["Absolute % error with confidence intervals (no smoothing)", "", "", ""];

plot_errors(
        times(3:length(times)), 
        errors(3:length(times)), 
        names(3:length(times)),
        titles; 
        stds(3:length(times)), 
        to_save_path="./visuals/January/diffusions/fixed_and_dynamic_diffusions_errors.png", 
        to_save=TO_SAVE
        )

work_precision_plot(names, algorithms, prob; DENSE = DENSE, 
                    SAVE_EVERYSTEP = false, to_save=TO_SAVE, 
                    to_save_path="./visuals/January/diffusions/fixed_diffusion_wp_diagram-comparing-the-orders.png",
                    title="EK1, fixed and dynamic diffusions comparison",
                    colors=[1 1 2 2 3 3 4 4 5 5]
                    )
