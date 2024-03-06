include("../../src/utils.jl")

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
evaltspan = (0.0, 50.0)
tstops = range(evaltspan..., length=5000)

# Generate the log-spaced values for rate 
exponents = range(-6, stop=5, length=150)
lengthscales = 10 .^ exponents

matern_prior_color = "red2"

algorithms = []
n_derivatives = 3
for l in lengthscales
        push!(algorithms, EK1(prior=Matern(n_derivatives, l), smooth=SMOOTH, diffusionmodel=DM))
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

solution_IWP = solve(prob, EK1(prior=IWP(n_derivatives), smooth=SMOOTH, diffusionmodel=DM), tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)
reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9, saveat=solution_IWP.t)

error_IWP = l2(solution_IWP.u, reference(solution_IWP.t))


# plotting the length scale vs Matern error
p = plot(layout=(1,1), legendfont = font(7), size = (800, 180))

xlabel = "length scale"

hline!( 
        p,
        [error_IWP],
        linestyle=:dash, linecolor=:black,
        label = "IWP",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [10mm 0mm],
        legend=:outerright,
        framestyle=:axes,
        fg_legend = :transparent
        )
plot!( 
        p,
        lengthscales[.!isnan.(all_errors)], 
        all_errors[.!isnan.(all_errors)],
        label = "Matern(length scale)",
        ylabel = "tRMSE", 
        xlabel = xlabel,
        title="",
        xaxis=:log10,
        left_margin = [20mm 0mm], 
        right_margin = [20mm 0mm],
        bottom_margin = [10mm 0mm],
        legend=:outerright,
        color = matern_prior_color,
        framestyle=:axes,
        fg_legend = :transparent,
        yaxis=:log10,
        dpi=600
        )

path = "./visuals/matern/length-scale-grid-search.pdf"
savefig(p, path)


# function find_first_min(array)
#     min_value = minimum(array[.!isnan.(array)])
#     return findfirst(x -> x == min_value, array)
# end

# print("Length scale:", find_first_min(all_errors))
# Length scale:139

# min_error = minimum(all_errors[.!isnan.(all_errors)])
# min_error_index = argmin(all_errors[.!isnan.(all_errors)])
# print("Minimum error of the Matern prior is $(min_error) at length scale $(lengthscales[min_error_index])\n")
# print("The error made by a solver with IWP prior is $(error_IWP)\n")

# Minimum error of the Matern prior is 0.11957508022607141 at length scale 366.23377139033727
# The error made by a solver with IWP prior is 0.11957161493488402