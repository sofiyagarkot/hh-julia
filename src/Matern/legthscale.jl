include("../utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

# Generate the log-spaced values for rate 
exponents = range(-0.5, stop=3, length=100)
lengthscales = 10 .^ exponents


algorithms = []
n_derivatives = 3
for l in lengthscales
        push!(algorithms, EK1(prior=Matern(n_derivatives, l), smooth=SMOOTH, diffusionmodel=DM))
end
push!(algorithms, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM))

all_errors_ch = Dict()
for i in 1:length(algorithms[1:end-1])
        algorithm =  algorithms[i]
        print("i $i --> prior $(algorithm.prior) \n")

        solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        solutions = [solution]
        errors, stds = absolute_errors(solutions, prob)
        for ch in 1:4
                errors_channel = sum(errors[1][ch])
                if !haskey(all_errors_ch, ch)
                        all_errors_ch[ch] = [errors_channel]
                else
                        push!(all_errors_ch[ch], errors_channel)
                end
        end
end

solution_IWP = solve(prob, algorithms[end], tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
solutions_IWP = [solution_IWP]
errors_IWP, stds_IWP = absolute_errors(solutions_IWP, prob)



p = plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))
titles = ["Log total absolute error, EK1(prior), dt=0.01, FixedDiffusion",         
        "", 
        "", 
        ""
        ]
# xlabels =["", "", "", "lengthscales"]
xlabels =["", "", "", "log10(lengthscales)"]
ylabels =["V", "m", "h", "n"]

for ch in 1:4
        plot!( 
                p[ch],
                log10.(lengthscales), 
                # lengthscales, 
                log.(all_errors_ch[ch]),
                label="Matern(3, l)",
                ylabel = ylabels[ch], 
                xlabel = xlabels[ch],
                title=titles[ch],
                left_margin = [20mm 0mm], 
                right_margin = [20mm 0mm],
                )
        log_sum_errors = log.(sum(errors_IWP[1][ch]))
        log_sum_errors_IWP = fill(log_sum_errors, length(lengthscales))

        plot!(
                p[ch],
                log10.(lengthscales),
                # lengthscales,
                log_sum_errors_IWP,
                label="IWP(3)",
                color=:black,
                ylabel = ylabels[ch], 
                xlabel = xlabels[ch],
                title=titles[ch],
                left_margin = [20mm 0mm], 
                right_margin = [20mm 0mm],
                )
end

display(p)


path = "./visuals/Matern/log-lengthscale-vs-total-abs-error.png"
savefig(p, path)
