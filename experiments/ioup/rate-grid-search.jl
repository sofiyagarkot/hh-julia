include("../../src/utils.jl")
SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

ioup_prior_color = "purple4"
# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

# exponents = range(0, stop=40, length=30)
# rates = -10 .^ exponents
# rates = vcat(rates, -rates)
# rates = sort(rates)
rates = range(-30, stop=30, length=50)

algorithms = []
n_derivatives = 3
for l in rates
        push!(algorithms, EK1(prior=IOUP(n_derivatives, l, update_rate_parameter=false), smooth=SMOOTH, diffusionmodel=DM))
end

all_errors = []
names_of_channels = ["V", "m", "h", "n"]

for i in 1:length(algorithms)
        algorithm =  algorithms[i]

        try
                solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
        catch e
                push!(all_errors, NaN)
                continue
        else
                solution = solve(prob, algorithm, tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)

                reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution.t)
                
                error = l2(solution.u, reference(solution.t))
                push!(all_errors, error)
        end
end

min_error = minimum(all_errors[.!isnan.(all_errors)])
min_error_index = argmin(all_errors[.!isnan.(all_errors)])

solution_IWP = solve(prob, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM), tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_IWP.t)

error_IWP = l2(solution_IWP.u, reference(solution_IWP.t))



p = plot(layout=(1,1), legendfont = font(7), size = (800, 180))

xlabel = "rate"

hline!( 
        p,
        [error_IWP],
        linestyle=:dash, linecolor=:black,
        label = "IWP",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        legend=:outerright,
        framestyle=:axes,
        fg_legend = :transparent
        )
plot!( 
        p,
        rates[.!isnan.(all_errors)], 
        all_errors[.!isnan.(all_errors)],
        label = "IOUP(rate)",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [5mm 0mm],
        legend=:outerright,
        color = "black",
        framestyle=:axes,
        fg_legend = :transparent,
        yaxis=:log10,
        )

scatter!(
        p, 
        [rates[min_error_index]],
        [min_error],
        label="minimum, rate=$(round(rates[min_error_index], digits=1))",
        left_margin = [20mm 0mm], 
        color = ioup_prior_color,
        right_margin = [20mm 0mm],
        mshape=:circle,
        markersize=3,
        yaxis=:log10,
        dpi=600

)

path = "./visuals/ioup/rate-grid-search.pdf"
savefig(p, path)

# print("Minimal error of IOUP(3, rate) at rate $(rates[min_error_index]) is $(min_error)\n")
# print("Error of IWP(3) is $error_IWP\n")
