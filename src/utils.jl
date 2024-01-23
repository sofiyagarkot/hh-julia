using ProbNumDiffEq, ComponentArrays
using UnPack
# using PlotlyJS, LaTeXStrings
using Plots, TuePlots
using OrdinaryDiffEq
using Statistics, LinearAlgebra
using Plots.PlotMeasures
using Colors
using DiffEqDevTools
using Distributions

Plots.theme(
    :dao;
    markerstrokewidth=1.2,
    legend=:left,
)
# Plots.theme(:default;
#     TuePlots.get_plotsjl_theme_kwargs(
#         TuePlots.SETTINGS[:AISTATS];
#         fontsize = true,
#         figsize = true,
#         single_column = true,
# )...)

function define_problem(duration = 50.0)
    αm(V, VT) = -0.32 * (V - VT - 13) / (exp(-(V - VT - 13) / 4) - 1)
    βm(V, VT) = 0.28 * (V - VT - 40) / (exp((V - VT - 40) / 5) - 1)

    αn(V, VT) = -0.032 * (V - VT - 15) / (exp(-(V - VT - 15) / 5) - 1)
    βn(V, VT) = 0.5 * exp(-(V - VT - 10) / 40)

    αh(V, VT) = 0.128 * exp(-(V - VT - 17) / 18)
    βh(V, VT) = 4 / (1 + exp(-(V - VT - 40) / 5))

    m_inf(V, VT) = 1 / (1 + βm(V, VT) / αm(V, VT))
    n_inf(V, VT) = 1 / (1 + βn(V, VT) / αn(V, VT))
    h_inf(V, VT) = 1 / (1 + βh(V, VT) / αh(V, VT))

    ENa = 53
    EK = -107
    area = 15e-5
    C = 1
    Eleak = -70
    VT = -60
    gleak = 0.1
    V0 = -70

    I(t) = (0 <= t <= duration) ? 500one(t) : zero(t)

    function f(du, u, params, t)
        V, m, n, h = u

        @unpack gNa, gK = params

        I_inj = 500*1e-6

        du[2] = dmdt = (αm(V, VT) * (1 - m) - βm(V, VT) * m)
        du[3] = dndt = (αn(V, VT) * (1 - n) - βn(V, VT) * n)
        du[4] = dhdt = (αh(V, VT) * (1 - h) - βh(V, VT) * h)

        INa = gNa * m^3 * h * (V - ENa) * area
        IK = gK * n^4 * (V - EK) * area
        Ileak = gleak * (V - Eleak) * area
        Cm = C * area
        du[1] = dVdt = -(Ileak + INa + IK - I_inj) / Cm
        # transform ODE to log the m,n  and h 
    end

    u0 = [V0, m_inf(V0, VT), n_inf(V0, VT), h_inf(V0, VT)]
    tspan = (0.0, duration)
    p = ComponentVector(gNa=20.0, gK=15.0)

    prob = ODEProblem(f, u0, tspan, p)
    return prob
end


function absolute_erros(solutions::Vector, prob::ODEProblem; baseline=Vern9())
    """
    Returns the absolute errors of the solutions a
    """
    all_errors = []
    all_stds = []
    
    for solution in solutions
        errors = []
        stds = []
        reference = solve(prob, baseline, abstol=1e-9, reltol=1e-9, saveat=solution.t)
        
        for ch in 1:4
            ref_ = reduce(hcat, reference.u)[ch, :]
            if isa(solution, ProbNumDiffEq.ProbODESolution)
                sol_ = reduce(hcat, mean.(solution.pu))[ch, :]
                push!(stds, reduce(hcat, std.(solution.pu))[ch, :])
                push!(errors, abs.(sol_ .- ref_))
            else
                sol_ = reduce(hcat, solution.u)[ch, :]
                push!(errors, abs.(sol_ .- ref_))
            end
        end
        
        push!(all_errors, errors)
        push!(all_stds, stds)
    end
    
    return all_errors, all_stds
end



function plot_errors_vs_dt(
    dts::Array, prob, algorithms, names; 
    to_save_path="./visuals/Januarydiffusions/errors-for-dts.png", 
    evaltspan = (0.0, 50.0), 
    TO_SAVE=false,
    colors=[],
    titles=["Average absolute % error, FixedDiffusion", "", "", ""],
    ylabels=["abs(error) (V), %", "abs(error) (m), %", "abs(error) (n), %", "abs(error) (h), %"],
    xlabels= ["", "", "", "step size, dt"],
    kwargs...
    )

    if isempty(colors)
        colors = 1:length(names)
    end

    p = plot(layout = (4, 1), 
            # legendfont = font(7), 
            size = (1000, 600)
            )
    for i in 1:length(algorithms)
        
        algorithm =  algorithms[i]
    
        for ch in 1:4
            alg_dts = [] 
            alg_errors = [] 
            for dt in dts
                n_samples = Int(maximum(evaltspan)/dt)
                tstops = range(evaltspan..., length=n_samples)
                try
                    solution = solve(prob, algorithm, tstops=tstops; kwargs...)
                    solutions = [solution]
                    errors = absolute_percentage_errors(solutions, prob)
                    errors = mean(errors[1][ch])
                    push!(alg_dts, dt)
                    push!(alg_errors, errors)
                catch PosDefException
                    print("\n\ti > ", i, " dt = ", dt)
                    continue
                end
            end

            plot!( 
                p[ch],
                log.(alg_dts), 
                log.(alg_errors), 
                label = names[i],
                legend=true,
                ylabel = ylabels[ch], 
                xlabel = xlabels[ch],
                title=titles[ch],
                color=colors[i],
                left_margin = [20mm 0mm], 
                right_margin = [20mm 0mm],
                legendfontsize=6,
                )
        end
        
    end


    if TO_SAVE
        savefig(p, to_save_path)
    else
        display(p)
    end

end

function plot_errors(
    times::Array, errors::Array, labels::Array, titles::Array; 
    stds = [], to_save_path="./visuals/sample.png", to_save=false, 
    colors = 1:length(labels),
    ylabels = ["abs error, V", "abs error, m", "abs error, n", "abs error, h"],
    xlabels = ["", "", "", "t"],
    legends = [false, false, false, true],
    p = nothing,
    log_plot = false,
    )
    if isa(p, Nothing)
        p = Plots.plot(layout = (4, 1), 
        # legendfont = font(7), 
                size = (1000, 600)
                )
    end

    for i in 1:length(labels)
        for j in 1:4

            if !isempty(stds)

                if log_plot
                    plot!(p[j], times[i], log.(errors[i][j]), xlabel = xlabels[j]; 
                            label = labels[i], color = colors[i], ylabel = ylabels[j], 
                            title=titles[j], legend = legends[j], left_margin = [20mm 0mm], 
                            right_margin = [20mm 0mm])

                    plot!(p[j], times[i], log.(stds[i][j]), linestyle=:dash, xlabel = xlabels[j]; 
                            label = labels[i]*" std", color = colors[i], ylabel = ylabels[j], 
                            title=titles[j], legend = legends[j], left_margin = [20mm 0mm], 
                            right_margin = [20mm 0mm])
                else
                    plot!(p[j], times[i], log.(errors[i][j]), xlabel = xlabels[j], ribbon = 3stds[i][j], 
                                fillalpha=0.2; label = labels[i], color = colors[i], ylabel = ylabels[j], 
                                title=titles[j], legend = legends[j], left_margin = [20mm 0mm], right_margin = [20mm 0mm])
                end
            else
                plot!(p[j], times[i], log.(errors[i][j]), xlabel = xlabels[j]; 
                            label = labels[i], color = colors[i], ylabel = ylabels[j], title=titles[j], 
                            legend = legends[j], left_margin = [20mm 0mm], right_margin = [20mm 0mm])
            end
        end
    end

    if to_save
        Plots.savefig(p, to_save_path)
    else
        display(p)
    end
    return p
end


function get_solutions(algorithms::Array, prob::ODEProblem; kwargs...)
    solutions = []
    for i in 1:length(algorithms)
        solution = solve(prob, algorithms[i]; kwargs...)
        push!(solutions, solution)
    end
    return solutions
end

function plot_solution_with_uncertainty(
    solutions::Array, derivative::Int, titles::Array, labels::Array;
     ts = [],
     log_std=false,
     to_save=false,
     to_save_path=""
     )
    ylabels = ["V", "m", "n", "h"]
    if isempty(ts)
        ts = solutions[1].t
    else
        ts = ts
    end

    # from Nathanael's banner
    interp(sol, ts) = StructArray([
        ProbNumDiffEq.interpolate(
            t,
            sol.t,
            sol.x_filt,
            sol.x_smooth,
            sol.diffusions,
            sol.cache;
            smoothed=sol.alg.smooth,
        )
        for t in ts
    ])
    p = plot(layout = (4, 1), legendfont = font(7), size = (1000, 600))

    for i in 1:length(solutions)
        solution = solutions[i]
        for ch in 1:4
            xs = interp(solution, ts)

            vecvec2mat(vv) = hcat(vv...)'
            
            if derivative == 0
                H = solution.cache.E0
            elseif derivative == 1
                H = solution.cache.E1
            elseif derivative == 2
                H = solution.cache.E2
            elseif derivative == 3
                H = solution.cache.E3
            end
                
            
            m = vecvec2mat(xs.μ) * H'
            std = sqrt.(vecvec2mat(diag.(xs.Σ))) * H'
            if log_std
                std_log = log.(std[:, ch])          
                ch_ylabel = ylabels[ch]      
                plot!(p[ch], ts, std_log, fillalpha = 0.5, 
                    xlabel = "t", ylabel = "log std "*"$ch_ylabel", 
                    title=titles[ch], legend=:topleft, label=labels[i],
                    left_margin = [20mm 0mm], right_margin = [20mm 0mm])

            else
                plot!(p[ch], ts, m[:, ch], ribbon = std[:, ch], 
                    fillalpha = 0.5, xlabel = "t", ylabel = "u(t)", 
                    title=titles[ch], legend=:topleft, label=labels[i], 
                    left_margin = [20mm 0mm], right_margin = [20mm 0mm])
            end
        end
    end
    if to_save
        savefig(p, to_save_path)
    else
        display(p)
    end
end

function work_precision_plot(
                            names::Array, 
                            algorithms :: Array, 
                            prob::ODEProblem; 
                            
                            abstols = [],
                            reltols = [],
                            DENSE = false, 
                            SAVE_EVERYSTEP = false, 
                            to_save=false, 
                            to_save_path="./visuals/sample_wp.png",
                            to_save_path2 = "./visuals/sample_wp2.png",
                            title = "Dynamic diffusion",
                            colors = [],
                            dts=[],
                            error_estimate = :l∞,
                            kwargs...
                            )
    _setups = []
    for i in 1:length(names)
        if !isempty(dts)
            push!(_setups, names[i] => Dict(:alg=>algorithms[i], :dts=>dts))
        else
            push!(_setups, names[i] => Dict(:alg=>algorithms[i]))
        end
    end

    labels = first.(_setups)
    setups = last.(_setups)

    test_sol = solve(prob, Vern9(), abstol=1/10^14, reltol=1/10^14)
    
    if isempty(colors)
        colors = 1:length(names)
    end

    if isempty(abstols)
        abstols = 1.0 ./ 10.0 .^ (6:10)
    end

    if isempty(reltols)
        reltols = 1.0 ./ 10.0 .^ (3:7)
    end

    wp = WorkPrecisionSet(
        prob, abstols, reltols, setups;
        names = labels,
        appxsol = test_sol,
        dense = DENSE,
        save_everystep = SAVE_EVERYSTEP,
        maxiters = Int(1e7),
        error_estimate = error_estimate,
        numruns = 5,
        kwargs...

    )

    p = plot(wp, title=title, color=colors, legend=:best, legendfontsize=8)

    if to_save
        savefig(p, to_save_path)
    else
        display(p)
    end

    p2 = plot(wp, x=:dts, y=:l∞, color=colors, title = title, legend=:best, legendfontsize=8)
    if to_save
        savefig(p2, to_save_path2)
    else
        display(p2)
    end
    return [p, p2]
end


function plot_2_work_precision_plots(
    names::Array, 
    algorithms :: Array, 
    names2::Array, 
    algorithms2 :: Array, 
    prob::ODEProblem; 
    
    abstols = [],
    reltols = [],
    DENSE = false, 
    SAVE_EVERYSTEP = false, 
    to_save=false, 
    to_save_path="./visuals/sample_wp.png",
    title = "Dynamic diffusion",
    dts = [],
    kwargs...
    )

    _setups = []
    for i in 1:length(names)
        if !isempty(dts)
            push!(_setups, names[i] => Dict(:alg=>algorithms[i], :dts=>dts))
        else
            push!(_setups, names[i] => Dict(:alg=>algorithms[i]))
        end
    end

    labels = first.(_setups)
    setups = last.(_setups)

    _setups2 = []
    for i in 1:length(names2)
        push!(_setups2, names2[i] => Dict(:alg=>algorithms2[i]))
    end

    labels2 = first.(_setups2)
    setups2 = last.(_setups2)


    test_sol = solve(prob, Vern9(), abstol=1/10^14, reltol=1/10^14)
    abstols = reltols = repeat([missing], length(dts))

    wp1 = WorkPrecisionSet(
        prob, abstols, reltols, setups;
        names = labels,
        appxsol = test_sol,
        dense = DENSE,
        save_everystep = SAVE_EVERYSTEP,
        maxiters = Int(1e7),
        numruns = 5,
        kwargs...
    )

    abstols = 1.0 ./ 10.0 .^ (6:10)
    reltols = 1.0 ./ 10.0 .^ (3:7)

    p = plot(wp1, title=title, legend=:topleft, legendfontsize=6)

    wp2 = WorkPrecisionSet(
        prob, abstols, reltols, setups2;
        names = labels2,
        appxsol = test_sol,
        dense = DENSE,
        save_everystep = SAVE_EVERYSTEP,
        maxiters = Int(1e7),
        numruns = 5,
    )

    plot!(p, wp2, x=:final, title=title, 
        fontsize=7, 
        legend=:topleft, 
        legendfontsize=6,
        )

    if to_save
        savefig(p, to_save_path)
    else
        display(p)
    end
end


function sample_prior(m0::Vector{Float64}, C0::Matrix{Float64}, ts::Vector{Float64}, Ah::Matrix{Float64} , Qh::Matrix{Float64}, d :: Int64, q::Int64 )
    sample = zeros(length(ts), d*(q+1))
    print(typeof(rand(MvNormal(m0, C0))), "\n", "shapr", size(rand(MvNormal(m0, C0))))
    sample[1, :] = rand(MvNormal(m0, C0))
    
    for i in 2:length(ts)
        sample[i, :] = rand(MvNormal(Ah*sample[i-1, :], Matrix(Qh)))
    end
    
    return sample
end

function plot_sample(
        ts::Vector{Float64}, sample::Matrix{Float64}, 
        E_all::Tuple, axes::Plots.Plot{Plots.GRBackend}, titles::Array, 
        derivative::Int64, color::Int64, dimension::Int64; label="")

    ylabels = ["V", "m", "n", "h"]

    E = E_all[derivative+1]
    Yi = sample * E'
    # color = parse(Colorant, "skyblue")
    for j in 1:dimension

        # ddt = derivative == 0 ? "" : L"$\frac{{d^{derivative}}}{{dt^{derivative}}}$"
        ddt = derivative == 0 ? "" : "d^$derivative/dt^$derivative"
        if isempty(label)
            
            if j == 4
                plot!(axes[j], ts, Yi[:, j], color=color, alpha = 0.5, label="",  legend=false, ylabel = "$ddt$(ylabels[j])", 
                    title=titles[j], xlabel="t", left_margin = [20mm 0mm], right_margin = [20mm 0mm])
            else 
                plot!(axes[j], ts, Yi[:, j], color=color, alpha=0.5, legend=false, ylabel = "$ddt$(ylabels[j])", 
                    title=titles[j],left_margin = [20mm 0mm], label="", right_margin = [20mm 0mm])
            end
        else
            if j == 4
                plot!(axes[j], ts, Yi[:, j], color=color, alpha = 0.5, label=label,legend=true, ylabel = "$ddt$(ylabels[j])", 
                    title=titles[j], xlabel="t", left_margin = [20mm 0mm], right_margin = [20mm 0mm])
            else 
                plot!(axes[j], ts, Yi[:, j], color=color,legend=false, alpha=0.5, ylabel = "$ddt$(ylabels[j])", 
                    title=titles[j],left_margin = [20mm 0mm], right_margin = [20mm 0mm])
            end
        end
    end
    return axes
end

