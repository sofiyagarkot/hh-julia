using TaylorSeries
include("utils.jl")

# defining inverse and forward transforms
identity =  x -> 1*x
deriv_logit = x -> ForwardDiff.derivative(logit, x)
deriv_identity = x -> ForwardDiff.derivative(identity, x)

forward_transforms = [identity, logit, logit, logit]
deriv_forward_transforms = [deriv_identity, deriv_logit, deriv_logit, deriv_logit]
inverse_transforms = [identity, sigmoid, sigmoid, sigmoid]

# plotting the transformed prior
tspan = range(0, 2, length=50)
prior = IWP(3)
p = plot(layout=(1,1), legendfont = font(7), size = (300, 180)) 
# plot(prior, tspan)
# color = "green3"

samples = ProbNumDiffEq.sample(prior, collect(tspan), 100)

all_0_derivative_samples = []
for t in samples
    zeros_derivative_samples = t[1, :]
    push!(all_0_derivative_samples, zeros_derivative_samples)
end

m = mean.(all_0_derivative_samples, dims=1)
std_ = std.(all_0_derivative_samples, dims=1)
m = vcat(m...)
std_ = vcat(std_...)

scale = 2
lower_bound = m - scale*std_
upper_bound = m + scale*std_

# plot!(p, tspan, inverse_transforms[2].(m), 
#         xlabel="time, s",label=L"$\mu$")
plot!(p, [0, 2], [0.5, 0.5], linecolor=:"green3", label=L"\mu", alpha=0.7)
       
plot!(p, tspan, 
        inverse_transforms[2].(lower_bound), 
        fillrange = inverse_transforms[2].(upper_bound), 
        fillalpha = 0.35,
        color="green3",
        alpha=0.25,
        xlabel="time, s",
        label=L"$\mu \pm$2σ",
        bottom_margin = [5mm 0mm],)

matrix_all_0_derivative_samples = hcat(all_0_derivative_samples...)
for s in eachrow(matrix_all_0_derivative_samples[1:15, :])
    plot!(p, tspan, inverse_transforms[2].(s), label="", color="green3", alpha = 0.3)
end
# hline!(p, [0.5], linecolor=:black, label=L"\mu", alpha=0.5)
plot!(p, 
        title="Transformed IWP(3) prior",
        legend=:bottomright,
        dpi=600)
p
savefig(p, "./visuals/transformed-prior.png")