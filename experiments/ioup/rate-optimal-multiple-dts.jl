include("../../src/utils.jl")
using  LaTeXStrings

# Colors
ioup_prior_color = "violet"
iwp_prior_color = "black"
ioup_optimal_color = "red"

dts = 10.0 .^ range(-1, -4, length=11)[1:end-1]
evaltspan = (0.0, 50.0)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()

rate_fixed_ioup = 12.9

prob = define_problem()

p = plot(legendfont = font(7), size = (500, 400))

# IWP(3) 
prior = IWP(3)
algorithm = EK1(prior=prior, smooth=DENSE, diffusionmodel=DM)
errors_iwp = []
for dt in dts
    try

        solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)

        reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
        error = l2(solution.u, reference(solution.t))
        push!(errors_iwp, error)

    catch e
        
        push!(errors_iwp, NaN)
        continue

    else            
    end

end

plot!(p, dts, errors_iwp, label="IWP", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_prior_color,
    fg_legend = :transparent)

scatter!(p, dts, errors_iwp, label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=iwp_prior_color,
    fg_legend = :transparent)
    


# IOUP(3) 

errors_ioup = []

for i in 1:length(dts)
    dt = dts[i]
    r = rate_fixed_ioup
    algorithm = EK1(prior=IOUP(3,  r, update_rate_parameter=false), smooth=DENSE, diffusionmodel=DM)
    try
        solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)

        reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
        error = l2(solution.u, reference(solution.t))
        push!(errors_ioup, error)
    catch e
        push!(errors_ioup, NaN)
        continue
    
    end

end

plot!(p, dts[.!isnan.(errors_ioup)], 
    float.(errors_ioup[.!isnan.(errors_ioup)]), 
    # linestyle=:dash, 
    label="IOUP($rate_fixed_ioup)", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=ioup_prior_color,
    fg_legend = :transparent)

scatter!(p, dts[.!isnan.(errors_ioup)], 
    float.(errors_ioup[.!isnan.(errors_ioup)]), label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=ioup_prior_color,
    fg_legend = :transparent)



# IOUPoptimal rate
function optimize_rate(rates, problem, dt)
    errors = []
    for rate in rates
        algorithm = EK1(prior=IOUP(3, rate, update_rate_parameter=false), smooth=DENSE, diffusionmodel=DM)
        try
            solution = solve(problem, algorithm, dt=dt, adaptive=false, dense = false)
            reference = solve(problem, Vern9(), abstol=1e-9, reltol=1e-9)
            error = l2(solution.u, reference(solution.t))
            push!(errors, error)
        catch e
            push!(errors, NaN)
            continue
        end
    end
    min_index = argmin(errors[.!isnan.(errors)])
    return rates[min_index]
end 

rates_optimal = []
errors_ioup = []
p
for i in 1:length(dts)
    dt = dts[i]
    rates = range(-300, 300, length=100)
    r = optimize_rate(rates, prob, dt)
    push!(rates_optimal, r)
    algorithm = EK1(prior=IOUP(3,  r, update_rate_parameter=false), smooth=DENSE, diffusionmodel=DM)
    try
        solution = solve(prob, algorithm, dt=dt, adaptive=false, dense = false)

        reference = solve(prob, Vern9(), abstol=1e-9, reltol=1e-9)
        error = l2(solution.u, reference(solution.t))
        push!(errors_ioup, error)
    catch e
        push!(errors_ioup, NaN)
        continue
    
    end

end

plot!(p, dts[.!isnan.(errors_ioup)], 
    float.(errors_ioup[.!isnan.(errors_ioup)]), 
    linestyle=:dash, 
    label="IOUP(optimal rate)", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=ioup_optimal_color,
    fg_legend = :transparent)

scatter!(p, dts[.!isnan.(errors_ioup)], 
    float.(errors_ioup[.!isnan.(errors_ioup)]), label="", 
    framestyle=:axes,
    xaxis=:log10, yaxis=:log10,
    color=ioup_optimal_color,
    fg_legend = :transparent)

plot!(p, xaxis=:log10, yaxis=:log10, 
        legend=:bottomright, xlabel=L"$\Delta$t", 
        ylabel="tRMSE",
        dpi=600)
p

savefig(p, "./visuals/ioup/rate-optimal-multiple-dts.pdf")