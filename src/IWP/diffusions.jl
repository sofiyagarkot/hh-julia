include("../utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false


prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

# Generate the log-spaced values for rate 
exponents = range(-100, stop=100, length=5)
initial_diffusions = 10 .^ exponents
algorithms = []

# fixed diffusion
DIFFUSION_MODELS = [
    FixedDiffusion(initial_diffusion = initial_diffusions[i], calibrate=false)
    for i in 1:length(initial_diffusions)
    ]
    
for diffusionmodel in DIFFUSION_MODELS
    push!(algorithms, EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=diffusionmodel))
end

all_errors_ch = Dict()
for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        solutions = [solution]
        errors, stds = absolute_errors(solutions, prob)
        for ch in 1:4
                errors_channel = sum(errors[1][ch])
                # errors_channel = mean(errors[1][ch])
                if !haskey(all_errors_ch, ch)
                        all_errors_ch[ch] = [errors_channel]
                else
                        push!(all_errors_ch[ch], errors_channel)
                end
        end
end



p = plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))

titles = ["Log total absolute error, EK1, dt=0.01, FixedDiffusion",         
        "", 
        "", 
        ""
        ]
xlabels =["", "", "", "log(diffusions)"]
ylabels =["V", "m", "h", "n"]

for ch in 1:4
        plot!( 
                p[ch],
                log10.(initial_diffusions), 
                log.(all_errors_ch[ch]),
                label="IWP(3)",
                ylabel = ylabels[ch], 
                xlabel = xlabels[ch],
                title=titles[ch],
                left_margin = [20mm 0mm], 
                right_margin = [20mm 0mm],
                )
        # scatter!(
        #         p[ch],  
        #         rates, 
        #         log.(all_errors_ch[ch]), 
        #         label="")
end

display(p)

# IOUP 

algorithms = []

# fixed diffusion
DIFFUSION_MODELS = [
    FixedDiffusion(initial_diffusion = initial_diffusions[i], calibrate=false)
    for i in 1:length(initial_diffusions)
    ]
    
for diffusionmodel in DIFFUSION_MODELS
    push!(algorithms, EK1(prior=IOUP(3, 14, update_rate_parameter=false), smooth=SMOOTH, diffusionmodel=diffusionmodel))
end

all_errors_ch = Dict()
for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        solutions = [solution]
        errors, stds = absolute_errors(solutions, prob)
        for ch in 1:4
                errors_channel = sum(errors[1][ch])
                # errors_channel = mean(errors[1][ch])
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
            log10.(initial_diffusions), 
            log.(all_errors_ch[ch]),
            label="IOUP(3, rate=14)",
            ylabel = ylabels[ch], 
            xlabel = xlabels[ch],
            title=titles[ch],
            left_margin = [20mm 0mm], 
            right_margin = [20mm 0mm],
            )
end

display(p)
# Matern

algorithms = []

# fixed diffusion
DIFFUSION_MODELS = [
    FixedDiffusion(initial_diffusion = initial_diffusions[i], calibrate=false)
    for i in 1:length(initial_diffusions)
    ]
    
for diffusionmodel in DIFFUSION_MODELS
    push!(algorithms, EK1(prior=Matern(3, 10), smooth=SMOOTH, diffusionmodel=diffusionmodel))
end

all_errors_ch = Dict()
for i in 1:length(algorithms)
        algorithm =  algorithms[i]

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

for ch in 1:4
    plot!( 
            p[ch],
            log10.(initial_diffusions), 
            log.(all_errors_ch[ch]),
            label="Matern(3, 10)",
            ylabel = ylabels[ch], 
            xlabel = xlabels[ch],
            title=titles[ch],
            left_margin = [20mm 0mm], 
            right_margin = [20mm 0mm],
            )
end



display(p)

path = "./visuals/fixed_diffusions/IWP-IOUP-Matern.png"
savefig(p, path)
