include("../../../src/utils.jl")
include("../../../src/customprior.jl")
using BlockDiagonals

SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = true

prob = define_problem()
rate_parameter = -1.0
num_derivatives = 3
dimension_iwp = 1
dimension_ioup = 3
dimension = dimension_iwp + dimension_ioup


F_iwp, L_iwp = get_drift_and_diffusion_matrices(IWP(dimension_iwp, num_derivatives))
F_ioup, L_ioup = get_drift_and_diffusion_matrices(IOUP(dimension_ioup, num_derivatives, rate_parameter))
# F_ioup, L_ioup = get_drift_and_diffusion_matrices(IWP(dimension_ioup, num_derivatives))
F_2, L_2 = get_drift_and_diffusion_matrices(IWP(dimension, num_derivatives))

F = Matrix(BlockDiagonal([F_iwp, F_ioup]))
L = Matrix(BlockDiagonal([L_iwp, L_ioup]))

# F, L = kron(I(2), [0.0 1.0; 0.0 0.0]), kron(I(2), [0.0; 1.0])
sol1 = solve(prob, EK1(order=1));
sol2 = solve(prob, EK1(prior=CustomSDEPrior(dimension, num_derivatives, F, L)));
names = ["IWP($num_derivatives)", "1xIWP($num_derivatives) + 3xIOUP($num_derivatives, $rate_parameter)"]  
solutions = [sol1, sol2] 


std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
titles=["Absolute % error with confidence intervals", "", "", ""];
plot_errors(times, errors, names, titles; stds, to_save_path="./visuals/January/different_priors/IWP_IOUP/errors_EK1.png", to_save=TO_SAVE)


plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/different_priors/IWP_IOUP/log_errors_EK1.png",
    to_save=false,
    ribbon = false
    )

errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
titles=["Absolute errors with confidence intervals", "", "", ""];
plot_errors(times, errors, names, titles; stds, to_save_path="./visuals/January/different_priors/IWP_IOUP/abs-errors_EK1.png", to_save=TO_SAVE)


plot_log_errors(
    times, errors, names, titles; 
    stds, 
    to_save_path="./visuals/January/different_priors/IWP_IOUP/log_abs_errors_EK1.png",
    to_save=TO_SAVE,
    ribbon = false
    )

