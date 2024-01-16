prob = define_problem()

sol_iwp3 = solve(prob, EK1(prior=IOUP(3, update_rate_parameter=true), smooth=true,), dt=0.01, adaptive=false)
sol_iwp4 = solve(prob, EK1(prior=IOUP(4, update_rate_parameter=true), smooth=true), dt=0.01, adaptive=false)
sol_iwp5 = solve(prob, EK1(prior=IOUP(5, update_rate_parameter=true), smooth=true), dt=0.01, adaptive=false)

solutions = [sol_iwp3, sol_iwp4, sol_iwp5]
times = [sol_iwp3.t, sol_iwp4.t, sol_iwp5.t]
std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
labels = ["IOUP3", "IOUP4", "IOUP5"]
stds = [stds_matern5, stds_ioup3, stds_ioup5, stds_iwp3, stds_iwp5]
titles=["Absolute % error with confidence intervals", "", "", ""];
plot_errors(times, errors, labels, titles; stds, to_save_path="./visuals/sample.png", to_save=false)
