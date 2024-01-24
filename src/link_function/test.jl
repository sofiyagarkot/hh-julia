include("../utils.jl")
using Plots

SMOOTH = DENSE = false
ADAPTIVE = false
TO_SAVE = false

# fixed diffusion
DM = FixedDiffusion()

prob = define_problem()
prob = define_problem_with_link_function()
evaltspan = (0.0, 1.0)
tstops = range(evaltspan..., length=10000)

algorithm = [EK1(prior=IWP(3), smooth=SMOOTH, diffusionmodel=DM)]
names = ["EK1(IWP(3))"]

solution = solve(prob, algorithms[1], tstops=tstops, adaptive=ADAPTIVE, dense=DENSE)

