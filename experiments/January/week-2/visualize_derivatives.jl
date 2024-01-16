prob = define_problem()
# EXPERIMENTING WITH HOW DOES SOLUTION LOOK LIKE 
SMOOTH = DENSE = false
names = [
        "EK1(IOUP3)",
        "EK1(Matern3)", 
        "EK1(IWP3)", 
        ]

algorithms = [
             EK1(prior=IOUP(3, update_rate_parameter=true), smooth=SMOOTH, diffusionmodel=FixedDiffusion()),
             EK1(prior=Matern(3, 1), smooth=SMOOTH, diffusionmodel=FixedDiffusion()), 
             EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=FixedDiffusion()), 
            #  EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=FixedDiffusion()), 
            #  EK1(prior=IWP(5), smooth=SMOOTH, diffusionmodel=FixedDiffusion())
             ]

prob = define_problem()

names = [
        "EK1(2)", 
        "EK1(3)", 
        # "EK1(4)"
        ]
algorithms = [
             EK1(prior=IWP(2), smooth=SMOOTH,diffusionmodel=FixedDiffusion()), 
             EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=FixedDiffusion()), 
        #      EK1(prior=IWP(4), smooth=SMOOTH, diffusionmodel=FixedDiffusion())
             ]
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=500)
solutions = get_solutions(algorithms, prob; tstops=tstops, dense = DENSE, adaptive=false)  

std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
titles=["Absolute % error with confidence intervals", "", "", ""];
plot_errors(times, errors, names, titles; stds, to_save_path="./visuals/sample.png", to_save=false)
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=500)
solutions = get_solutions(algorithms, prob; tstops=tstops, adaptive=false, dense=false)  

std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
titles=["Absolute % error with confidence intervals", "", "", ""];
plot_errors(times, errors, names, titles; stds, to_save_path="./visuals/January/sample.png", to_save=false)

SMOOTH = false
algorithms = [
             EK1(prior=IOUP(3, update_rate_parameter=true), smooth=SMOOTH, ),
             EK1(prior=Matern(3, 1), smooth=SMOOTH), 
             EK1(prior=IWP(3), smooth=SMOOTH, ), 
            #  EK1(prior=IWP(4), smooth=true, diffusionmodel=FixedDiffusion()), 
            #  EK1(prior=IWP(5), smooth=true, diffusionmodel=FixedDiffusion())
             ]
work_precision_plot(names, algorithms, prob; DENSE = false, SAVE_EVERYSTEP = false, to_save=true, to_save_path="./visuals/January/wp_diagram-comparing-the-orders.png")


ch = 1
evaltspan = (0.0, 55.0)
ts = range(evaltspan..., length=300)
plot_solution_with_uncertainty( solutions[1], ch, 0, "IOUP: Voltage", ts=ts)
plot_solution_with_uncertainty( solutions[2], ch, 0, "Matern 3,1: Voltage", ts=ts)
plot_solution_with_uncertainty( solutions[3], ch, 0, "IWP: Voltage", ts=ts)

plot_solution_with_uncertainty( solutions[1], ch, 1, "IOUP: Voltage", ts=ts)
plot_solution_with_uncertainty( solutions[2], ch, 1, "Matern 3,1: Voltage", ts=ts)
plot_solution_with_uncertainty( solutions[3], ch, 1, "IWP: Voltage", ts=ts)

plot_solution_with_uncertainty( solutions[1], ch, 2, "IOUP: Voltage", ts=ts)
plot_solution_with_uncertainty( solutions[2], ch, 2, "Matern 3,1: Voltage", ts=ts)
plot_solution_with_uncertainty( solutions[3], ch, 2, "IWP: Voltage", ts=ts)


# MORE
algorithms = [
             EK1(prior=IOUP(3, update_rate_parameter=true), smooth=true, diffusionmodel=FixedDiffusion()),
             EK1(prior=Matern(3, 1), smooth=true, diffusionmodel=FixedDiffusion()), 
             EK1(prior=Matern(3, 2), smooth=true, diffusionmodel=FixedDiffusion()), 
             EK1(prior=Matern(3, 3), smooth=true, diffusionmodel=FixedDiffusion()), 
             EK1(prior=IWP(3), smooth=true, diffusionmodel=FixedDiffusion()), 
            #  EK1(prior=IWP(4), smooth=true, diffusionmodel=FixedDiffusion()), 
            #  EK1(prior=IWP(5), smooth=true, diffusionmodel=FixedDiffusion())
             ]