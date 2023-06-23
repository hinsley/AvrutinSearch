module AvrutinSearch
export homoclinic_to_equilibrium

using IntervalArithmetic
using Peaks
using Roots


function homoclinic_to_equilibrium(f, x_star::Float64, T_0::Interval{Float64}, I::Interval{Float64}, r_max::UInt64, n_iter::UInt64=0, tol::Float64=1e-5)::Union{Bool, Nothing}
    # Return values:
    # true: Connection from x_star to T found.
    # false: Connection from x_star to T does not exist.
    # nothing: Connection from x_star to T not found.

    # Check if x_star is a fixed point.
    if f(x_star) != x_star
        throw(ArgumentError("x_star is not a fixed point."))
    end

    # Check if x_star is in T_0 and I.
    if !(x_star ∈ T_0) || !(x_star ∈ I)
        throw(ArgumentError("x_star must be in T_0 and I."))
    end
    
    # Check if f(I) = I.
    if f(I) != I
        throw(ArgumentError("f(I) must equal I."))
    end

    # Compute T as `n_iter'th iterate of T_0.
    T = T_0
    for _ in 1:n_iter
        T = f(T)
    end

    # Compute k, the intervals V_1, ..., V_k and functions f_1^{-1}, ..., f_k^{-1}.
    sample_points = I.lo:tol:I.hi
    sample_values = f.(sample_points)
    partition = sort(vcat(Peaks.argmaxima(sample_values), Peaks.argminima(sample_values)))
    V = [Interval(sample_points[partition[i]], sample_points[partition[i+1]]) for i in 1:2:length(partition)-1]
    k = length(V)
    f_inv = [y -> find_zero((x -> f(x) - y), (V[i].lo, V[i].hi)) for i in 1:k]

    # Proceed with the algorithm from the paper.
    S_p = [(x_star, 0)]
    i_max = 0
    while length(S_p) > 0 # While S_p is not empty.
        p, i = popfirst!(S_p)
        if i <= r_max
            i_max = max(i_max, i)
            for j in 1:k
                if p ∈ V[j]
                    x = f_inv[j](p)
                    if x ∈ T && x != x_star
                        return true
                    end
                    if x ∈ I && x != x_star
                        push!(S_p, (x, i+1))
                    end
                end
            end
        end
    end

    if i_max == r_max
        return false
    else
        return nothing
    end
end

end # module AvrutinSearch