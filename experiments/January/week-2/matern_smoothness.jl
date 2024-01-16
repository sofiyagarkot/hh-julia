include("../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = true

prob = define_problem()

names = [
        "Matern(1, 1)",
        "Matern(2, 1)", 
        "Matern(3, 1)", 
        "Matern(4, 1)", 
        "Matern(5, 1)"
        ]
SMOOTH = false

# Comparing the dimensions of an ODE

algorithms = [
             EK1(prior=Matern(1, 1), smooth=SMOOTH),
             EK1(prior=Matern(2, 1), smooth=SMOOTH), 
             EK1(prior=Matern(3, 1), smooth=SMOOTH), 
             EK1(prior=Matern(4, 1), smooth=SMOOTH), 
             EK1(prior=Matern(5, 1), smooth=SMOOTH)
             ]

# DM = FixedDiffusion()
# algorithms = [
#         EK1(prior=Matern(1, 1), smooth=SMOOTH, diffusionmodel=DM),
#         EK1(prior=Matern(2, 1), smooth=SMOOTH, diffusionmodel=DM), 
#         EK1(prior=Matern(3, 1), smooth=SMOOTH, diffusionmodel=DM), 
#         EK1(prior=Matern(4, 1), smooth=SMOOTH, diffusionmodel=DM), 
#         EK1(prior=Matern(5, 1), smooth=SMOOTH, diffusionmodel=DM)
#         ]

evaltspan = (0.0, 50.0)
# tstops = range(evaltspan..., length=5000)
solutions = get_solutions(algorithms, prob; adaptive=ADAPTIVE, dense=DENSE)  
titles=["Error with Matern prior (lengthscale = 1, no smoothing, non-adaptive steps, dynamic diffusion)", "", "", ""];

std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]

plot_errors(
        times, errors, names, titles; 
        stds, to_save_path="./visuals/January/matern/abs-error-different-orders.png", 
        to_save=TO_SAVE
        )

plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/matern/log-abs-error-different-orders.png", 
        to_save=TO_SAVE,
        ribbon = false
        )

names_3_5 = [ 
        "Matern(3, 1)", 
        "Matern(4, 1)", 
        "Matern(5, 1)"
        ]
SMOOTH = false

# Comparing the dimensions of an ODE

algorithms_3_5 = [
             EK1(prior=Matern(3, 1), smooth=SMOOTH), 
             EK1(prior=Matern(4, 1), smooth=SMOOTH), 
             EK1(prior=Matern(5, 1), smooth=SMOOTH)
             ]

solutions_3_5 = get_solutions(algorithms_3_5, prob; adaptive=ADAPTIVE, dense=DENSE)  
titles=["Error with Matern prior (lengthscale = 1, no smoothing, non-adaptive steps, dynamic diffusion)", "", "", ""];
times_3_5 =[sol_i.t for sol_i in solutions_3_5]
std_bool=true
# errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
errors_3_5, stds_3_5 = absolute_errors(solutions_3_5, prob, std_bool = std_bool)
plot_log_errors(
        times_3_5, errors_3_5, names_3_5, titles; 
        stds=stds_3_5, 
        to_save_path="./visuals/January/matern/log-abs-error-order-3-5.png", 
        to_save=TO_SAVE,
        ribbon = false
        )

# abstols = 1.0 ./ 10.0 .^ (3:8)
# reltols = 1.0 ./ 10.0 .^ (1:5)

# work_precision_plot(
#         names, algorithms, prob; DENSE = DENSE,
#         abstols=abstols, reltols=reltols,
#         SAVE_EVERYSTEP = false, to_save=TO_SAVE, 
#         to_save_path="./visuals/January/matern/wp_diagram-comparing-the-orders.png",
#         title = "Comparison of different orders of Matern prior (lengthscale=1)")

# Comparison of lengthscale parameters
dim = 3
names = [
        "Matern($dim, 1)",
        "Matern($dim, 2)", 
        "Matern($dim, 3)", 
        "Matern($dim, 4)", 
        "Matern($dim, 5)", 
        "Matern($dim, 6)", 
        "Matern($dim, 7)", 
        "Matern($dim, 8)",
        ]
SMOOTH = false

algorithms = [
             EK1(prior=Matern(dim, 1), smooth=SMOOTH, ),
             EK1(prior=Matern(dim, 2), smooth=SMOOTH, ), 
             EK1(prior=Matern(dim, 3), smooth=SMOOTH, ), 
             EK1(prior=Matern(dim, 4), smooth=SMOOTH, ), 
             EK1(prior=Matern(dim, 5), smooth=SMOOTH, ), 
             EK1(prior=Matern(dim, 6), smooth=SMOOTH, ), 
             EK1(prior=Matern(dim, 7), smooth=SMOOTH, ), 
             EK1(prior=Matern(dim, 8), smooth=SMOOTH, ), 
             ]

evaltspan = (0.0, 50.0)
solutions = get_solutions(algorithms, prob;  adaptive=ADAPTIVE, dense=DENSE)  
titles=["Error with Matern prior (dimension = $dim, no smoothing, adaptive steps, dynamic diffusion)", "", "", ""];

std_bool=true
errors, stds = absolute_percentage_errors(solutions, prob, std_bool = std_bool)
# errors, stds = absolute_errors(solutions, prob, std_bool = std_bool)
times = [sol_i.t for sol_i in solutions]
titles=["Absolute % error with Matern prior (dimension = $dim, no smoothing, adaptive steps, dynamic diffusion)", "", "", ""];

plot_errors(
        times, errors, names, titles; 
        stds, to_save_path="./visuals/January/matern/percentage-error-different-lengthscales.png", 
        to_save=TO_SAVE
        )
titles=["Log absolute % error with Matern prior (dimension = $dim, no smoothing, adaptive steps, dynamic diffusion)", "", "", ""];

plot_log_errors(
        times, errors, names, titles; 
        stds, 
        to_save_path="./visuals/January/matern/log-abs-percentage-error-different-lengthscales-ribbon.png", 
        to_save=TO_SAVE,
        ribbon = true
        )

# work_precision_plot(names, algorithms, prob; DENSE = DENSE, 
#                     SAVE_EVERYSTEP = false, to_save=TO_SAVE, 
#                     to_save_path="./visuals/January/matern/wp_diagram-comparing-the-lengthscales.png",
#                     title = "Comparison of different lengthscales of Matern prior (dimension=$dim)")

titles = ["Log mean absolute error",         
          "", 
          "", 
          ""
          ]
xlabels =["", "", "", "lengthscale"]
ylabels =["V", "m", "h", "n"]


lengthscales = 1:8
p = plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))
for ch in 1:4
        errors_ch = []
        for i in 1:length(algorithms)
                algorithm =  algorithms[i]
                solution = solve(prob, algorithm, dense=false)
                solutions = [solution]
                errors = absolute_errors(solutions, prob)
                errors = mean(errors[1][ch])
                push!(errors_ch, errors)
        end

        plot!( 
                p[ch],
                lengthscales, 
                log.(errors_ch),
                legend=false,
                ylabel = ylabels[ch], 
                xlabel = xlabels[ch],
                title=titles[ch],
                left_margin = [20mm 0mm], 
                right_margin = [20mm 0mm],
                legendfontsize=6,
                )
end
TO_SAVE = true
to_save_path = "./visuals/January/matern/log-abs-error-different-lengthscales-per-channel.png"
if TO_SAVE
        savefig(p, to_save_path)
else
        display(p)
end