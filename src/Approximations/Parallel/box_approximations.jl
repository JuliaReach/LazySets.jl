export box_approximation_symmetric_parallel

function box_approximation_symmetric_parallel(S::LazySet{N}
                                    )::Hyperrectangle{N} where {N<:Real}
    (c, r) = box_approximation_helper_parallel(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

@inline function box_approximation_helper_parallel(S::LazySet{N}) where {N<:Real}
    n = dim(S)
    c = SharedVector{N}(n)
    r = SharedVector{N}(n)

    distribute_task!(c, r, S)
    return convert(Array,c), convert(Array,r)
end

# Here's the kernel
function process_chunk!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}, irange::UnitRange{Int64}) where {N<:Real}

    d = zeros(N, dim(S))

    for i in irange
        d[i] = one(N)
        htop = ρ(d, S)
        d[i] = -one(N)
        hbottom = -ρ(d, S)
        d[i] = zero(N)
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
    end
end

"""
    symmetric_interval_hull_parallel

Alias for `box_approximation_symmetric_parallel`.
"""
symmetric_interval_hull_parallel = box_approximation_symmetric_parallel