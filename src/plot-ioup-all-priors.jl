include("utils.jl")

priors = [
    # IWP(3),
    IOUP(3, -30.0, update_rate_parameter=false),
    IOUP(3, 12.9, update_rate_parameter=false),
    IOUP(3, 20.0, update_rate_parameter=false),
    ]
labels = [
    # "IWP",
    "IOUP(-30)",
    "IOUP(12.9)",
    "IOUP(20)",
]
    
plotrange = range(0, 1.0, length=1000)

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
    # if i == 3
    #     color = "purple4"
    # else
    #     color = "black"
    # end

    plot!(p[i], priors[i], plotrange; title=labels[i], 
                # ylims=(-20, 20), 
                ylabel="", color="black",  dpi=600)
end
p

xlabel!(p[1], "time, ms")
xlabel!(p[2], "time, ms")
xlabel!(p[3], "time, ms")

# savefig(p, "./visuals/highlighted_ioup_priors.png")
savefig(p, "./visuals/essay/ioup_priors.png")
