using Pkg
Pkg.activate("AvrutinSearch")

using Interpolations
using IntervalArithmetic
using Roots

using AvrutinSearch


xs = collect(0.0:1e-7:1.0)
f(x) = 3.8 * (x - x^2)
fs = f.(xs)
f_interpolated = linear_interpolation(xs, fs)

I = 0..1
x_star = find_zero(x -> f_interpolated(x) - x, (0.5, 1))
T_0 = x_star ± 1e-1
r_max = UInt64(10)
n_iter = UInt64(0)

x, i = AvrutinSearch.homoclinic_to_equilibrium(xs, fs, x_star, T_0, I, r_max, n_iter)
println("x = $x, i = $i")

using Plots

indices = round.(Int, LinRange(1, length(xs), 100))
plt = plot(Shape([T_0.lo, T_0.hi, T_0.hi, T_0.lo], [I.lo, I.lo, I.hi, I.hi]), opacity=0.2, label="\$\\mathcal{T}_0\$", xlims=(I.lo, I.hi), ylims=(I.lo, I.hi), aspect_ratio=:equal)
plot!(plt, [I.lo, I.hi], [I.lo, I.hi], label="Bisectrix")
plot!(plt, xs[indices], fs[indices], label="f")
x_forward_iterates = [x]
for j in 1:i
    push!(x_forward_iterates, f_interpolated(x_forward_iterates[end]))
end
# Make points for cobweb plot of forward iterates.
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
plot!(plt, x_cobweb, y_cobweb, label="Forward iterates")
display(plt)