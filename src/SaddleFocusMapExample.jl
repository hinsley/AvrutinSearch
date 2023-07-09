using Pkg
Pkg.activate("AvrutinSearch")

using Interpolations
using IntervalArithmetic
using Plots
using Roots

using AvrutinSearch

@gif for μ in LinRange(-0.05, 0.27, 1000)
    xs = collect(-1.0:1e-5:1.0)
    ρ = 0.9
    ω = 10.0
    lw = 5.0 # Line width.
    f(x) = x == 0.0 ? 0.0 : sign(x) * (μ + abs(x)^ρ*cos(ω*log(abs(x))))
    fs = f.(xs)
    f_interpolated = linear_interpolation(xs, fs)

    I = -0.9..0.9
    x_star = find_zeros(x -> f_interpolated(x) - x, -1e-3, 1e-3)[1]
    T_0 = x_star ± 3e-4
    r_max = UInt64(3)
    n_iter = UInt64(0)

    #indices = round.(Int, LinRange(1, length(xs), 10000))
    #plt = plot(xs[indices], fs[indices], label="f", xlims=(I.lo, I.hi), ylims=(I.lo, I.hi), aspect_ratio=:equal)
    #xlims!(plt, (-0.95, 0.95))
    #ylims!(plt, (-0.95, 0.95))

    @time result = AvrutinSearch.homoclinic_to_equilibrium(xs, fs, x_star, T_0, I, r_max, n_iter)
    if result === nothing
        println("Connection from x_star to T not found.")

        indices = round.(Int, LinRange(1, length(xs), 10000))
        # Put legend at top left.
        plt = plot(Shape([T_0.lo, T_0.hi, T_0.hi, T_0.lo], [I.lo, I.lo, I.hi, I.hi]), opacity=0.2, label="\$\\mathcal{T}_0\$", title="\$\\mu=$(μ)\$", legend=:topleft, xlims=(I.lo, I.hi), ylims=(I.lo, I.hi), aspect_ratio=:equal, size=(800, 800))
        # Compute T as `n_iter'th iterate of T_0.
        T = T_0
        for i in 1:n_iter
            T = AvrutinSearch.iterate_interval(T, xs, fs, f)
            plot!(plt, Shape([T.lo, T.hi, T.hi, T.lo], [I.lo, I.lo, I.hi, I.hi]), opacity=0.2, label="\$\\mathcal{T}_$i\$")
        end
        plot!(plt, [I.lo, I.hi], [I.lo, I.hi], lw=lw, label="Bisectrix")
        plot!(plt, xs[indices], fs[indices], lw=lw, label="f")
        # Plot cobweb plot of forward iterates.
        plot!(plt, [], [], lw=lw, label="Forward iterates")
        display(plt)
    else
        x, i = result
        println("x = $x, i = $i")

        indices = round.(Int, LinRange(1, length(xs), 10000))
        plt = plot(Shape([T_0.lo, T_0.hi, T_0.hi, T_0.lo], [I.lo, I.lo, I.hi, I.hi]), opacity=0.2, label="\$\\mathcal{T}_0\$", title="\$\\mu=$(μ)\$", legend=:topleft, xlims=(I.lo, I.hi), ylims=(I.lo, I.hi), aspect_ratio=:equal, size=(800, 800))
        # Compute T as `n_iter'th iterate of T_0.
        T = T_0
        for i in 1:n_iter
            T = AvrutinSearch.iterate_interval(T, xs, fs, f)
            plot!(plt, Shape([T.lo, T.hi, T.hi, T.lo], [I.lo, I.lo, I.hi, I.hi]), opacity=0.2, label="\$\\mathcal{T}_$i\$")
        end
        plot!(plt, [I.lo, I.hi], [I.lo, I.hi], lw=lw, label="Bisectrix")
        plot!(plt, xs[indices], fs[indices], lw=lw, label="f")
        # Make points for cobweb plot of forward iterates.
        x_forward_iterates = [x]
        for j in 1:i
            push!(x_forward_iterates, f_interpolated(x_forward_iterates[end]))
        end
        x_cobweb = [x_forward_iterates[1]]
        y_cobweb = [x_forward_iterates[1]]
        for j in 1:i
            push!(x_cobweb, x_forward_iterates[j])
            push!(y_cobweb, x_forward_iterates[j])
            push!(x_cobweb, x_forward_iterates[j])
            push!(y_cobweb, x_forward_iterates[j+1])
        end
        push!(x_cobweb, x_star)
        push!(y_cobweb, x_star)
        # Plot cobweb plot of forward iterates.
        plot!(plt, x_cobweb, y_cobweb, lw=lw, label="Forward iterates")
        display(plt)
    end
end