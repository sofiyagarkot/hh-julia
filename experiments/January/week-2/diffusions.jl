# trying to see what does fixed diffusion do 
prob = define_problem()

evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

sol_ioup3_update = solve(prob, EK1(prior=IWP(3), smooth=true, diffusionmodel=FixedDiffusion()), tstops=tstops, adaptive=false)
sol_ioup3_dynamic = solve(prob, EK1(prior=IWP(3), smooth=true, diffusionmodel=DynamicDiffusion()), tstops=tstops, adaptive=false)
# sol_ioup3_dont_update = solve(prob, EK1(prior=IWP(3), smooth=true, diffusionmodel=DynamicDiffusion()), dt=0.01, adaptive=false)

solutions = [sol_ioup3_update, sol_ioup3_dynamic]
times = [sol_ioup3_update.t, sol_ioup3_dynamic.t]
std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
labels = ["FixedDiffusion", "DynamicDiffusion"]
titles=["Absolute % error with confidence intervals", "", "", ""];
plot_errors(times, errors, labels, titles; stds, to_save_path="./visuals/January/diffusions.png", to_save=false)