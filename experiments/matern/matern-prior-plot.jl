include("../../src/utils.jl")

priors = [
    Matern(3, 0.1),
    Matern(3, 1),
    Matern(3, 10)
    ]
labels = [
    "Matérn(0.1)",
    "Matérn(1)",
    "Matérn(10)"
]
    
plotrange = range(0, 1.0, length=1500)

p = plot(
    layout=(1, 3),
    legendfont = font(7),
    size = (800, 200),
    bottom_margin = [7mm 0mm],
    left_margin = [5mm 0mm],
    right_margin = [5mm 0mm],
    top_margin = [5mm 0mm],
    dpi=600
)

for i in 1:length(priors)

    color = "black"

    plot!(p[i], priors[i], plotrange; title=labels[i], 
                ylabel="", color=color,  dpi=600)
end

ylims!(p[1], (-3.1*(10^(-5)), 3.1*(10^(-5))))
ylims!(p[2], (-0.09, 0.09))
ylims!(p[3], (-270, 270))

xlabel!(p[1], "time, ms")
xlabel!(p[2], "time, ms")
xlabel!(p[3], "time, ms")

savefig(p, "./visuals/matern/matern-prior-plot.pdf")