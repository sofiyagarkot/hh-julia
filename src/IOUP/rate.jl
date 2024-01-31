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
rates = range(-30, stop=50, length=100)

algorithms = []
n_derivatives = 3
for l in rates
        push!(algorithms, EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false), smooth=SMOOTH, diffusionmodel=DM))
end

all_errors_ch = Dict()
all_errors_ch_spans = Dict()
names_of_channels = ["V", "m", "h", "n"]

for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        solutions = [solution]
        errors, stds = absolute_errors(solutions, prob)
        for ch in 1:4
                name_ch = names_of_channels[ch]
                errors_channel = sum(errors[1][ch])
                
                if !haskey(all_errors_ch, name_ch)
                        all_errors_ch[name_ch] = [errors_channel]
                else
                        push!(all_errors_ch[name_ch], errors_channel)
                end
        end
end

for name_ch in names_of_channels
        min_ch = minimum(all_errors_ch[name_ch])
        max_ch = maximum(all_errors_ch[name_ch])
        all_errors_ch_spans[name_ch] = log.([min_ch, max_ch])
end

push!(algorithms, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM))
solution_IWP = solve(prob, algorithms[end], tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
solutions_IWP = [solution_IWP]
errors_IWP, stds_IWP = absolute_errors(solutions_IWP, prob)
log_sum_errors = log.(sum(errors_IWP[1][1]))
log_sum_errors_IWP = fill(log_sum_errors, length(rates))

p = plot(layout=(1,1), legendfont = font(7), size = (800, 180))

min_v = round(all_errors_ch_spans["V"][1], digits=1)
max_v =round( all_errors_ch_spans["V"][2], digits=1)
min_m = round(all_errors_ch_spans["m"][1], digits=1)
max_m = round(all_errors_ch_spans["m"][2], digits=1)
min_h = round(all_errors_ch_spans["h"][1], digits=1)
max_h = round(all_errors_ch_spans["h"][2], digits=1)
min_n = round(all_errors_ch_spans["n"][1], digits=1)
max_n = round(all_errors_ch_spans["n"][2], digits=1)

long_label = "Error scales:\nV: $min_v - $max_v;\nm: $min_m - $max_m;\nh: $min_h - $max_h;\nn: $min_n - $max_n"

xlabel = "rate"

hline!( 
        p,
        [log_sum_errors],
        linestyle=:dash, linecolor=:black,
        label = "IWP(3)",
        ylabel = "", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        yticks=[],
        framestyle=:axes,
        fg_legend = :transparent
        )
plot!( 
        p,
        rates, 
        log.(all_errors_ch["V"]),
        label = "IOUP(3, rate)",
        ylabel = "", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        yticks=[],
        framestyle=:axes,
        fg_legend = :transparent
        )


p = annotate!(30, 6.1, text(long_label, 7, :center))
display(p)
path = "./visuals/IOUP/rate-vs-total-abs-error-neg-pos.png"
savefig(p, path)


# Log-plot

# Generate the log-spaced values for rate 
exponents = range(-1, stop=1.5, length=100)
rates = 10 .^ exponents
rates = sort(rates)


algorithms = []
n_derivatives = 3
for l in rates
        push!(algorithms, EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false), smooth=SMOOTH, diffusionmodel=DM))
end
# push!(algorithms, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM))

all_errors_ch = Dict()
all_errors_ch_spans = Dict()
names_of_channels = ["V", "m", "h", "n"]

for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        solutions = [solution]
        errors, stds = absolute_errors(solutions, prob)
        for ch in 1:4
                name_ch = names_of_channels[ch]
                errors_channel = sum(errors[1][ch])
                
                if !haskey(all_errors_ch, name_ch)
                        all_errors_ch[name_ch] = [errors_channel]
                else
                        push!(all_errors_ch[name_ch], errors_channel)
                end
        end
end

for name_ch in names_of_channels
        min_ch = minimum(all_errors_ch[name_ch])
        max_ch = maximum(all_errors_ch[name_ch])
        all_errors_ch_spans[name_ch] = log.([min_ch, max_ch])
end

push!(algorithms, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM))
solution_IWP = solve(prob, algorithms[end], tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
solutions_IWP = [solution_IWP]
errors_IWP, stds_IWP = absolute_errors(solutions_IWP, prob)
log_sum_errors = log.(sum(errors_IWP[1][1]))
log_sum_errors_IWP = fill(log_sum_errors, length(rates))

p = plot(layout=(1,1), legendfont = font(7), size = (800, 180))
# p=annotate!(0.5, 0.5, text(long_label, 7, :center))

min_v = round(all_errors_ch_spans["V"][1], digits=1)
max_v =round( all_errors_ch_spans["V"][2], digits=1)
min_m = round(all_errors_ch_spans["m"][1], digits=1)
max_m = round(all_errors_ch_spans["m"][2], digits=1)
min_h = round(all_errors_ch_spans["h"][1], digits=1)
max_h = round(all_errors_ch_spans["h"][2], digits=1)
min_n = round(all_errors_ch_spans["n"][1], digits=1)
max_n = round(all_errors_ch_spans["n"][2], digits=1)

long_text = "Error scales:\nV: $min_v - $max_v;\nm: $min_m - $max_m;\nh: $min_h - $max_h;\nn: $min_n - $max_n"

xlabel = "rate"

hline!( 
        p,
        [log_sum_errors],
        linestyle=:dash, linecolor=:black,
        label = "IWP(3)",
        ylabel = "", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        yticks=[],
        framestyle=:axes,
        fg_legend = :transparent
        )
plot!( 
        p,
        rates, 
        log.(all_errors_ch["V"]),
        label = "IOUP(3, rate)",
        ylabel = "", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        yticks=[],
        framestyle=:axes,
        fg_legend = :transparent
        )
# Put text outside of a plot
# Add text to the plot

minimums = [0.8686868686868687, 1.121212121212121]
indices = findall(x -> x in minimums, log10.(rates))
     
ms = [10^0.8686868686868687, 10^1.121212121212121]

scatter!(
        p,
        ms, 
        log.(all_errors_ch["V"])[indices],
        label="minimum",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        mshape=:circle,
        color=:orange,
        markersize=2
        )

p = annotate!(30, 6.1, text(long_label, 7, :center))

display(p)

path = "./visuals/IOUP/rate-vs-total-abs-error-pos.png"
savefig(p, path)


