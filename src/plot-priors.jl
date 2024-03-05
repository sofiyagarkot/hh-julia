include("utils.jl")

priors = [
    IWP(3),
    IOUP(3, -1000, update_rate_parameter=false),
    IOUP(3, 1.0, update_rate_parameter=false),
    Matern(3, 0.01),
    Matern(3, 1),
    Matern(3, 100)
    ]
labels = [
    "IWP(3)",
    "IOUP(3, -1000)",
    "IOUP(3, 1)",
    "Matern(3, 0.01)",
    "Matern(3, 1)",
    "Matern(3, 100)"
]
    
plotrange = range(0, 10, length=250)

p = plot(
    layout=(2, 3),
    legendfont = font(7),
    size = (800, 400),
    bottom_margin = [5mm 0mm],
    left_margin = [5mm 0mm],
    right_margin = [5mm 0mm],
    top_margin = [5mm 0mm],
    dpi=600
)

for i in 1:length(priors)
    plot!(p[i], priors[i], plotrange; title=labels[i], 
                ylabel="", color="black",  dpi=600)
end
p
savefig(p, "./visuals/all_priors.png")