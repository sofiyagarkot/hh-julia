include("utils.jl")

prob = define_problem()
min_dt = 0.001

DM = FixedDiffusion()
integ = init(prob, EK1(prior=IWP(3), diffusionmodel = DM), dt=min_dt, adaptive=false)

mean_all_V = []
mean_all_m = []

std_all_V = []
std_all_m = []

diffusions = []

E0 = integ.cache.E0 # just extract it from the integrator to make sure the sizes etc are correct
x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
initial_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
initial_stddev = std(x_0) 

integ.cache.local_diffusion 
integ.cache.global_diffusion
integ.cache.default_diffusion
default_diffusion = integ.cache.default_diffusion
push!(diffusions, default_diffusion)

push!(mean_all_V, initial_mean[1])
push!(std_all_V, initial_stddev[1])


push!(mean_all_m, initial_mean[2])
push!(std_all_m, initial_stddev[2])

tspan =  collect(0.0:min_dt:0.01)

for i in tspan[2:end]

    step!(integ)

    x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
    step_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
    step_stddev = std(x_0) 

    push!(mean_all_V, step_mean[1])
    push!(std_all_V, step_stddev[1])

    push!(mean_all_m, step_mean[2])
    push!(std_all_m, step_stddev[2])

    global_diffusion = integ.cache.global_diffusion
    push!(diffusions, global_diffusion)

end


non_transformed_m = mean_all_m
non_transformed_V = mean_all_V

non_transformed_std_m = diffusions .* std_all_m
non_transformed_std_V = diffusions .* std_all_V

p = plot(layout=(2,2), 
    legendfont = font(7), size = (700, 400)
    )
# plot!(p[1], tspan,  mean_all_V, ribbon = (10^14)*(diffusions .* std_all_V), label="V, ribbon scaled by 10^13", color="blue", fillalpha=0.2)

plot!(p[1], tspan,  mean_all_V, label="mean V", color="blue", fillalpha=0.2)

# plot!(p[3], tspan,  std_all_V, label="std V", color="blue", fillalpha=0.2)
plot!(p[3], tspan,  diffusions .* std_all_V, label="std V, \nscaled by diffusion", color="blue", fillalpha=0.2)
# plot!(p[2], tspan,  mean_all_m, ribbon = (10^13)*(diffusions .* std_all_m), label="m", color="blue", fillalpha=0.2)

plot!(p[2], tspan,  mean_all_m, label="mean m", color="blue", fillalpha=0.2)
# plot!(p[4], tspan,  std_all_m, label="std m", color="blue", fillalpha=0.2)
plot!(p[4], tspan,  diffusions .* std_all_V, label="std m, \nscaled by diffusion", color="blue", fillalpha=0.2)

title!(p[1], dpi = 600, left_margin = [10mm 0mm],  "Non-transformed space: mean and std of 1-step interpolation", titleloc = :left)
xlabel!(p[3], "time, ms", bottom_margin = [10mm 0mm],)
xlabel!(p[4], "time, ms")
ylabel!(p[1], "V, mV", left_margin = [10mm 0mm])
ylabel!(p[2], "probability")


path1 = "./visuals/extrapolation/extrapolation-2/no-ribbon-non-transformed.png"
savefig(p, path1)



prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V()
integ = init(prob_transformed, EK1(prior=IWP(3), diffusionmodel = DM), dt=min_dt, adaptive=false)

mean_all_V = []
mean_all_m = []

std_all_V = []
std_all_m = []

diffusions = []

E0 = integ.cache.E0 # just extract it from the integrator to make sure the sizes etc are correct
x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
initial_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
initial_stddev = std(x_0) 

integ.cache.local_diffusion 
integ.cache.global_diffusion
integ.cache.default_diffusion
default_diffusion = integ.cache.default_diffusion
push!(diffusions, default_diffusion)

push!(mean_all_V, initial_mean[1])
push!(std_all_V, initial_stddev[1])

push!(mean_all_m, initial_mean[2])
push!(std_all_m, initial_stddev[2])

tspan =  collect(0.0:min_dt:0.01)

for i in tspan[2:end]

    step!(integ)

    x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
    step_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
    step_stddev = std(x_0) 

    push!(mean_all_V, step_mean[1])
    push!(std_all_V, step_stddev[1])

    push!(mean_all_m, step_mean[2])
    push!(std_all_m, step_stddev[2])

    global_diffusion = integ.cache.global_diffusion
    push!(diffusions, global_diffusion)
end


p = plot(layout=(2,2), 
    legendfont = font(7), size = (700, 400)
    )
# plot!(p[1], tspan,  mean_all_V, ribbon = (10^14)*(diffusions .* std_all_V), label="V, ribbon scaled by 10^13", color="blue", fillalpha=0.2)

plot!(p[1], tspan,  mean_all_V, label="mean V", color="blue", fillalpha=0.2)

# plot!(p[3], tspan,  std_all_V, label="std V", color="blue", fillalpha=0.2)
plot!(p[3], tspan,  diffusions .* std_all_V, label="std V, \nscaled by diffusion", color="blue", fillalpha=0.2)


# plot!(p[2], tspan,  mean_all_m, ribbon = (10^13)*(diffusions .* std_all_m), label="m", color="blue", fillalpha=0.2)
plot!(p[2], tspan,  mean_all_m, label="mean m", color="blue", fillalpha=0.2)
# plot!(p[4], tspan,  std_all_m, label="std m", color="blue", fillalpha=0.2)
plot!(p[4], tspan,  diffusions .* std_all_V, label="std m, \nscaled by diffusion", color="blue", fillalpha=0.2)

title!(p[1], dpi = 600,"Transformed space: mean and std of 1-step interpolation", titleloc = :left)
xlabel!(p[3], "time, ms", bottom_margin = [10mm 0mm],)
xlabel!(p[4], "time, ms")
ylabel!(p[1], "", left_margin = [10mm 0mm])
ylabel!(p[2], "")

path2 = "./visuals/extrapolation/extrapolation-2/no-ribbon-transformed.png"
savefig(p, path2)

# p = plot(layout=(3,1), legendfont = font(7), size = (600, 400))
# plot!(p[1], tspan,  inverse_transforms[1].(mean_all_V), ribbon = inverse_transforms[1].((10^13)*(diffusions .* std_all_V)), label="V", color="blue", fillalpha=0.2, title="Mean and std of prediction in inversely transformed space")
# plot!(p[2], tspan,  inverse_transforms[1].(std_all_V), label="std V", color="blue", fillalpha=0.2)
# plot!(p[3], tspan,  inverse_transforms[1].(diffusions .* std_all_V), label="diffusion scaled std V", color="blue", fillalpha=0.2)


p = plot(layout=(2,2), 
    legendfont = font(7), size = (700, 400)
    )
# plot!(p[1], tspan,  inverse_transforms[1].(mean_all_V), ribbon = inverse_transforms[1].((10^14)*(diffusions .* std_all_V)), label="V, ribbon scaled by 10^13", color="blue", fillalpha=0.2)
plot!(p[1], tspan,  inverse_transforms[1].(mean_all_V), label="mean V", color="blue", fillalpha=0.2)
# plot!(p[3], tspan,  std_all_V, label="std V", color="blue", fillalpha=0.2)
plot!(p[3], tspan,  inverse_transforms[1].(diffusions .* std_all_V), label="std V, \nscaled by diffusion", color="blue", fillalpha=0.2)


# plot!(p[2], tspan,  inverse_transforms[2].(mean_all_m), ribbon = inverse_transforms[1].((10^13)*(diffusions .* std_all_m)), label="m", color="blue", fillalpha=0.2)
plot!(p[2], tspan,  inverse_transforms[2].(mean_all_m), label="mean m", color="blue", fillalpha=0.2)
# plot!(p[4], tspan,  std_all_m, label="std m", color="blue", fillalpha=0.2)
plot!(p[4], tspan,  inverse_transforms[2].(diffusions .* std_all_V), label="std m, \nscaled by diffusion", color="blue", fillalpha=0.2)

title!(p[1],dpi = 600, "Inversely transformed space: mean and std of 1-step interpolation", titleloc = :left)
xlabel!(p[3], "time, ms", bottom_margin = [10mm 0mm],)
xlabel!(p[4], "time, ms")
ylabel!(p[1], "V, mV", left_margin = [10mm 0mm])
ylabel!(p[2], "probability")

path3 = "./visuals/extrapolation/extrapolation-2/no-ribbon-inverse-transformed.png"
savefig(p, path3)

transformed_m = inverse_transforms[2].(mean_all_m)
transformed_V = inverse_transforms[1].(mean_all_V)

transformed_std_m = inverse_transforms[2].(diffusions .* std_all_m)
transformed_std_V = inverse_transforms[1].(diffusions .* std_all_V)



p = plot(layout=(2,2), 
    legendfont = font(7), size = (1200, 400)
    )
# plot!(p[1], tspan,  inverse_transforms[1].(mean_all_V), ribbon = inverse_transforms[1].((10^14)*(diffusions .* std_all_V)), label="V, ribbon scaled by 10^13", color="blue", fillalpha=0.2)
plot!(p[1], tspan,  transformed_V - non_transformed_V, label="inv_transformed(mean(V)) \n- mean(V)", color="blue", legend=:outerleft, fillalpha=0.2)
# plot!(p[3], tspan,  std_all_V, label="std V", color="blue", fillalpha=0.2)
plot!(p[3], tspan,  transformed_std_V - non_transformed_std_V, label="inv_transformed(diffusion*std(V))\n  - diffusion*std(V)", color="blue", legend=:outerleft, fillalpha=0.2)


# plot!(p[2], tspan,  inverse_transforms[2].(mean_all_m), ribbon = inverse_transforms[1].((10^13)*(diffusions .* std_all_m)), label="m", color="blue", fillalpha=0.2)
plot!(p[2], tspan,  transformed_m - non_transformed_m, label="inv_transformed(mean(m)) \n- mean(m)", legend=:outerright,color="blue", fillalpha=0.2)
# plot!(p[4], tspan,  std_all_m, label="std m", color="blue", fillalpha=0.2)
plot!(p[4], tspan,  transformed_std_m - non_transformed_std_m, label="inv_transformed(diffusion*std(m))\n  - diffusion*std(m)", color="blue", legend=:outerright, fillalpha=0.2)

title!(p[1],dpi = 600, "Difference between inversely transformed and non-transformed", titleloc = :left)
xlabel!(p[3], "time, ms", bottom_margin = [10mm 0mm],)
xlabel!(p[4], "time, ms")
# ylabel!(p[1], "V, mV", )
# ylabel!(p[2], "probability", )

path4 = "./visuals/extrapolation/extrapolation-2/no-ribbon-difference.png"
savefig(p, path4)
p


freq1 = 0.001
sin_x = x -> 0.5*sin(freq1*x) + 0.5
arcsin_x = x-> (1/freq1)*asin(2*x - 1)
deriv_arcsin =  x -> ForwardDiff.derivative(arcsin_x, x)

freq = 0.0001
arcsin_V = x -> (1/freq)*asin( (x+20)/90)
sin_v = x -> 90*sin(freq*x) - 20
deriv_arcsin_V =  x -> ForwardDiff.derivative(arcsin_V, x)

forward_transforms = [arcsin_V, arcsin_x, arcsin_x, arcsin_x]
deriv_forward_transforms = [deriv_arcsin_V, deriv_arcsin, deriv_arcsin, deriv_arcsin]
inverse_transforms = [sin_v, sin_x, sin_x, sin_x]

prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V(forward_transforms, deriv_forward_transforms, inverse_transforms)
integ = init(prob_transformed, EK1(prior=IWP(3), diffusionmodel = DM), dt=min_dt, adaptive=false)

mean_all_V = []
mean_all_m = []

std_all_V = []
std_all_m = []

diffusions = []

E0 = integ.cache.E0 # just extract it from the integrator to make sure the sizes etc are correct
x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
initial_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
initial_stddev = std(x_0) 

integ.cache.local_diffusion 
integ.cache.global_diffusion
integ.cache.default_diffusion
default_diffusion = integ.cache.default_diffusion
push!(diffusions, default_diffusion)

push!(mean_all_V, initial_mean[1])
push!(std_all_V, initial_stddev[1])

push!(mean_all_m, initial_mean[2])
push!(std_all_m, initial_stddev[2])

tspan =  collect(0.0:min_dt:0.01)

for i in tspan[2:end]

    step!(integ)

    x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
    step_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
    step_stddev = std(x_0) 

    push!(mean_all_V, step_mean[1])
    push!(std_all_V, step_stddev[1])

    push!(mean_all_m, step_mean[2])
    push!(std_all_m, step_stddev[2])

    global_diffusion = integ.cache.global_diffusion
    push!(diffusions, global_diffusion)
end

transformed_m_sin = inverse_transforms[2].(mean_all_m)
transformed_V_sin = inverse_transforms[1].(mean_all_V)

transformed_std_m_sin = inverse_transforms[2].(diffusions .* std_all_m)
transformed_std_V_sin = inverse_transforms[1].(diffusions .* std_all_V)



# p = plot(layout=(2,2), 
#     legendfont = font(7), size = (1200, 400)
#     )
# plot!(p[1], tspan,  inverse_transforms[1].(mean_all_V), ribbon = inverse_transforms[1].((10^14)*(diffusions .* std_all_V)), label="V, ribbon scaled by 10^13", color="blue", fillalpha=0.2)
plot!(p[1], tspan,  transformed_V_sin - non_transformed_V, label="sin inv_transformed(mean(V)) \n- mean(V)", color="red", legend=:outerleft, fillalpha=0.2)
# plot!(p[3], tspan,  std_all_V, label="std V", color="blue", fillalpha=0.2)
plot!(p[3], tspan,  transformed_std_V_sin - non_transformed_std_V, label="sin inv_transformed(diffusion*std(V))\n  - diffusion*std(V)", color="red", legend=:outerleft, fillalpha=0.2)


# plot!(p[2], tspan,  inverse_transforms[2].(mean_all_m), ribbon = inverse_transforms[1].((10^13)*(diffusions .* std_all_m)), label="m", color="blue", fillalpha=0.2)
plot!(p[2], tspan,  transformed_m_sin - non_transformed_m, label="sin inv_transformed(mean(m)) \n- mean(m)", legend=:outerright,color="red", fillalpha=0.2)
# plot!(p[4], tspan,  std_all_m, label="std m", color="blue", fillalpha=0.2)
plot!(p[4], tspan,  transformed_std_m_sin - non_transformed_std_m, label="sin inv_transformed(diffusion*std(m))\n  - diffusion*std(m)", color="red", legend=:outerright, fillalpha=0.2)

title!(p[1],dpi = 600, "Difference between inversely transformed and non-transformed", titleloc = :left)
xlabel!(p[3], "time, ms", bottom_margin = [10mm 0mm],)
xlabel!(p[4], "time, ms")
# ylabel!(p[1], "V, mV", )
# ylabel!(p[2], "probability", )

path4 = "./visuals/extrapolation/extrapolation-2/no-ribbon-sin-difference.png"
savefig(p, path4)




# on sinusoids with frequencies 1.0 and 1.0 it performs the same as IWP(3)
# with bigger frequencies -- worse

freq1 = 5.0
sin_x = x -> 0.5*sin(freq1*x) + 0.5
arcsin_x = x-> (1/freq1)*asin(2*x - 1)
deriv_arcsin =  x -> ForwardDiff.derivative(arcsin_x, x)

freq = 5.0
arcsin_V = x -> (1/freq)*asin( (x+20)/90)
sin_v = x -> 90*sin(freq*x) - 20
deriv_arcsin_V =  x -> ForwardDiff.derivative(arcsin_V, x)

forward_transforms = [arcsin_V, arcsin_x, arcsin_x, arcsin_x]
deriv_forward_transforms = [deriv_arcsin_V, deriv_arcsin, deriv_arcsin, deriv_arcsin]
inverse_transforms = [sin_v, sin_x, sin_x, sin_x]

prob_transformed, forward_transforms, deriv_forward_transforms, inverse_transforms, ode_func = generate_transformed_settings_V(forward_transforms, deriv_forward_transforms, inverse_transforms)
integ = init(prob_transformed, EK1(prior=IWP(3), diffusionmodel = DM), dt=min_dt, adaptive=false)

mean_all_V = []
mean_all_m = []

std_all_V = []
std_all_m = []

diffusions = []

E0 = integ.cache.E0 # just extract it from the integrator to make sure the sizes etc are correct
x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
initial_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
initial_stddev = std(x_0) 

integ.cache.local_diffusion 
integ.cache.global_diffusion
integ.cache.default_diffusion
default_diffusion = integ.cache.default_diffusion
push!(diffusions, default_diffusion)

push!(mean_all_V, initial_mean[1])
push!(std_all_V, initial_stddev[1])

push!(mean_all_m, initial_mean[2])
push!(std_all_m, initial_stddev[2])

tspan =  collect(0.0:min_dt:0.01)

for i in tspan[2:end]

    step!(integ)

    x_0 = E0 * integ.cache.x # this gives you a Gaussian that describes the zero'th derivative of x
    step_mean = mean(x_0) # this computes the standard deviation of the Gaussian; equivalent to sqrt.(diag(x_0.Σ))
    step_stddev = std(x_0) 

    push!(mean_all_V, step_mean[1])
    push!(std_all_V, step_stddev[1])

    push!(mean_all_m, step_mean[2])
    push!(std_all_m, step_stddev[2])

    global_diffusion = integ.cache.global_diffusion
    push!(diffusions, global_diffusion)
end

transformed_m_sin = inverse_transforms[2].(mean_all_m)
transformed_V_sin = inverse_transforms[1].(mean_all_V)

transformed_std_m_sin = inverse_transforms[2].(diffusions .* std_all_m)
transformed_std_V_sin = inverse_transforms[1].(diffusions .* std_all_V)



# p = plot(layout=(2,2), 
#     legendfont = font(7), size = (1200, 400)
#     )
# plot!(p[1], tspan,  inverse_transforms[1].(mean_all_V), ribbon = inverse_transforms[1].((10^14)*(diffusions .* std_all_V)), label="V, ribbon scaled by 10^13", color="blue", fillalpha=0.2)
plot!(p[1], tspan,  transformed_V_sin - non_transformed_V, label="sin 5.0 inv_transformed(mean(V)) \n- mean(V)", color="purple", legend=:outerleft, fillalpha=0.2)
# plot!(p[3], tspan,  std_all_V, label="std V", color="blue", fillalpha=0.2)
plot!(p[3], tspan,  transformed_std_V_sin - non_transformed_std_V, label="sin 5.0 inv_transformed(diffusion*std(V))\n  - diffusion*std(V)", color="purple", legend=:outerleft, fillalpha=0.2)


# plot!(p[2], tspan,  inverse_transforms[2].(mean_all_m), ribbon = inverse_transforms[1].((10^13)*(diffusions .* std_all_m)), label="m", color="blue", fillalpha=0.2)
plot!(p[2], tspan,  transformed_m_sin - non_transformed_m, label="sin 5.0 inv_transformed(mean(m)) \n- mean(m)", legend=:outerright,color="purple", fillalpha=0.2)
# plot!(p[4], tspan,  std_all_m, label="std m", color="blue", fillalpha=0.2)
plot!(p[4], tspan,  transformed_std_m_sin - non_transformed_std_m, label="sin 5.0 inv_transformed(diffusion*std(m))\n  - diffusion*std(m)", color="purple", legend=:outerright, fillalpha=0.2)

title!(p[1],dpi = 600, "Difference between inversely transformed and non-transformed", titleloc = :left)
xlabel!(p[3], "time, ms", bottom_margin = [10mm 0mm],)
xlabel!(p[4], "time, ms")
# ylabel!(p[1], "V, mV", )
# ylabel!(p[2], "probability", )