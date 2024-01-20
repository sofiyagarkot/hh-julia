include("../utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = true
TO_SAVE = false

# TODO: make a scatter plot for the narrower grid 


prob = define_problem()

exponents = range(0, stop=2, length=20)

# Generate the log-spaced values
rates = 10 .^ exponents
rates = vcat(rates, -rates)
rates = sort(rates)
# rates = range(-40, stop=-10, length=2)
p = plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))
for n_derivatives in 1:5
    algorithms = []
    for l in rates
            push!(algorithms, EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false), smooth=SMOOTH))
    end


    titles = ["Log mean absolute error, EK1(IOUP(update=false))",         
            "", 
            "", 
            ""
            ]
    xlabels =["", "", "", "rates"]
    ylabels =["V", "m", "h", "n"]

    all_errors_ch = Dict()
    for i in 1:length(algorithms)
            algorithm =  algorithms[i]
            print("i $i --> prior $(algorithm.prior) \n")

            solution = solve(prob, algorithm, dense=false)
            solutions = [solution]
            errors = absolute_errors(solutions, prob)
            for ch in 1:4
                    errors_channel = mean(errors[1][ch])
                    if !haskey(all_errors_ch, ch)
                            all_errors_ch[ch] = [errors_channel]
                    else
                            push!(all_errors_ch[ch], errors_channel)
                    end
            end
    end
    for ch in 1:4
            plot!( 
                    p[ch],
                    rates, 
                    log.(all_errors_ch[ch]),
                    label="IOUP($n_derivatives)",
                    ylabel = ylabels[ch], 
                    xlabel = xlabels[ch],
                    title=titles[ch],
                    left_margin = [20mm 0mm], 
                    right_margin = [20mm 0mm],
                    legendfontsize=6,
                    )
    end
end
display(p)

TO_SAVE = true
to_save_path = "./visuals/January/ioup/log-abs-error-different-rates-and-n-derivatives-100-100.png"
if TO_SAVE
        savefig(p, to_save_path)
else
        display(p)
end