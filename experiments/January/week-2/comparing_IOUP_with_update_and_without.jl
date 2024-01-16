include("../../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = false


prob = define_problem()

names = [
        "IOUP(3), update=true",
        "IOUP(3, -1), update=false",
        "IOUP(3, 1), update=false",
        "IWP(3)"
        ]
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=500)
algorithms = [
        EK1(prior=IOUP(3,  update_rate_parameter=true)),
        EK1(prior=IOUP(3, -1, update_rate_parameter=false)),
        EK1(prior=IOUP(3, 1, update_rate_parameter=false)),
            #  EK1(prior=IOUP(3, update_rate_parameter=true), smooth=SMOOTH),
            #  EK1(prior=IOUP(3, -1, update_rate_parameter=false), smooth=SMOOTH),
        EK1(prior=IWP(3), smooth=SMOOTH),
        ]


# gettiong solutions for adaptive steps and dynamic diffusion
solutions = get_solutions(
    algorithms, prob; 
    dense = DENSE, adaptive=true,
    tstops=tstops
    )  

    
std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
# errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
titles=["Absolute % error EK1 (dimension = 3, no smoothing, adaptive steps, dynamic diffusion)", "", "", ""];

plot_errors(
        times, errors, names, titles; 
        stds, to_save_path="./visuals/January/ioup/percentage-errors.png", 
        to_save=true
        )

titles=["Log absolute % error EK1 (dimension = 3, no smoothing, adaptive steps, dynamic diffusion)", "", "", ""];

plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/ioup/log-abs-percentage-error-ribbon.png", 
        to_save=true,
        ribbon = true
        )
        
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
titles=["Log absolute error EK1 (dimension = 3, no smoothing, adaptive steps, dynamic diffusion)", "", "", ""];

plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/ioup/log-abs-error.png", 
    to_save=true,
    ribbon = false
    )


ch = 1
evaltspan = (0.0, 56.1)
ts = range(evaltspan..., length=5000)
derivative = 1
titles = ["Solutions for"*" $derivative derivative", "", "", ""]

plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        log_std=false,
        to_save=true,
        to_save_path="./visuals/January/ioup/solutions-fixed_diffusion_1st_voltage.png"
        )

titles = ["Log std for"*" $derivative derivative", "", "", ""]

plot_solution_with_uncertainty(
    solutions, derivative, titles, names;
        ts = ts,
        log_std=true,
        to_save=true,
        to_save_path="./visuals/January/ioup/log_std-fixed_diffusion_0th_voltage.png"
        )