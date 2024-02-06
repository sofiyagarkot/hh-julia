include("../utils.jl")

dts = 10.0 .^ range(-2, -4, length=5)

DENSE = SMOOTH = TO_SAVE = false
ADAPTIVE = false
SAVE_EVERYSTEP = false

DM = FixedDiffusion()
evaltspan = (0.0, 50.0)

_setups = [
    # "Implicit Euler" => Dict(:alg=>ImplicitEuler(), :dts=>dts)
    "Exponential Euler" => Dict(:alg=>LawsonEuler(krylov = true, m = 50), :dts=>dts) # on krylov and m : https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/ 
    "EK1(1)" => Dict(:alg=>EK1(order=1, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    "EK1(2)" => Dict(:alg=>EK1(order=2, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    "EK1(3)" => Dict(:alg=>EK1(order=3, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    "EK1(4)" => Dict(:alg=>EK1(order=4, smooth=DENSE, diffusionmodel=DM), :dts=>dts)
    ]
    
labels = first.(_setups)
setups = last.(_setups)

# COLORS
colors = [colorant"darkorange"]
# push!(colors, colorant"purple")

c1 = colorant"lightblue"
c2 = colorant"darkblue"
colors2 = range(c1, stop=c2, length=4)

for c in colors2
    push!(colors, c)
end


prob = define_problem()

test_sol = solve(prob, Vern9(), abstol=1/10^14, reltol=1/10^14)
abstols = reltols = repeat([missing], length(dts))

wp = WorkPrecisionSet(
        prob, abstols, reltols, setups;
        names = labels,
        appxsol = test_sol,
        dense = DENSE,
        save_everystep = SAVE_EVERYSTEP,
        maxiters = Int(1e7),
        numruns = 5,
        error_estimate = :l2,
        adaptive = ADAPTIVE,
    )



p = plot(wp, x=:dts,y=:l2, title="", 
    fontsize=7, 
    legend=:bottomright, 
    color=colors',
    legendfontsize=8,
    )
p
# savefig(p, "./visuals/baseline/fixed_diffusion_wp_EK1_IWP_ExpEuler.png")
