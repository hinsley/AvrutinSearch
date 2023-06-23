module AvrutinSearch
export homoclinic_to_equilibrium

using Interpolations
using IntervalArithmetic
using Peaks
using Roots


function homoclinic_to_equilibrium(xs::Vector{Float64}, fs::Vector{Float64}, x_star::Float64, T_0::Interval{Float64}, I::Interval{Float64}, r_max::UInt64, n_iter::UInt64=0)::Union{Tuple{Float64, UInt64}, Nothing}
    # Return values:
    # true: Connection from x_star to T found.
    # false: Connection from x_star to T does not exist.
    # nothing: Connection from x_star to T not found.

    # Interpolate f from f_pairs.
    f = linear_interpolation(xs, fs)

    # Check if x_star is a fixed point.
    if f(x_star) != x_star
        throw(ArgumentError("x_star is not a fixed point."))
    end

    # Check if x_star is in T_0 and I.
    if !(x_star ∈ T_0)
        throw(ArgumentError("x_star must be in T_0."))
    end

    # Check if T_0 is a subset of I.
    if !(T_0 ⊆ I)
        throw(ArgumentError("T_0 must be a subset of I."))
    end
    
    # TODO: Implement iterates of intervals for linearly interpolated functions.
    # Check if f(I) = I.
    #if f(I) != I
    #    throw(ArgumentError("f(I) must equal I."))
    #end

    # Compute T as `n_iter'th iterate of T_0.
    T = T_0
    #for _ in 1:n_iter
    #    T = f(T)
    #end

    # Compute k, the intervals V_1, ..., V_k and functions f_1^{-1}, ..., f_k^{-1}.
    sample_points = vcat([I.lo], filter(x -> I.lo < x < I.hi, xs), [I.hi])
    sample_values = f.(sample_points)
    partition = sort(vcat(1, Peaks.argmaxima(sample_values), Peaks.argminima(sample_values), length(sample_points)))
    V = [Interval(sample_points[partition[i]], sample_points[partition[i+1]]) for i in 1:length(partition)-1]
    k = length(V)
    # Get the ranges of f over each V_i.
    f_V = []
    for i in 1:k
        push!(f_V, Interval(sort([f(V[i].lo), f(V[i].hi)])...))
    end
    f_inv = []
    for i in 1:k
        # Produce a linearly interpolated local inverse of f.
        # We have to supply the domain points in increasing order for Interpolations.jl.
        if sample_values[partition[i]] < sample_values[partition[i+1]] # f is increasing on V_i.
            f_inv_i = linear_interpolation(sample_values[partition[i]:partition[i+1]], sample_points[partition[i]:partition[i+1]])
        else # f is decreasing on V_i.
            f_inv_i = linear_interpolation(sample_values[partition[i+1]:-1:partition[i]], sample_points[partition[i+1]:-1:partition[i]])
        end
        push!(f_inv, f_inv_i)
    end

    # Proceed with the algorithm from the paper.
    S_p = [(x_star, UInt64(0))]
    i_max = UInt64(0)
    while length(S_p) > 0 # While S_p is not empty.
        p, i = popfirst!(S_p)
        if i <= r_max
            i_max = max(i_max, i)
            for j in 1:k
                if p ∈ f_V[j]
                    x = f_inv[j](p)
                    if x ∈ T && x != x_star
                        return (x, i+1)
                    end
                    if x ∈ I && x != x_star
                        push!(S_p, (x, i+1))
                    end
                end
            end
        end
    end

    if i_max == r_max
        return nothing
    else
        return (NaN, UInt64(0))
    end
end

end # module AvrutinSearch