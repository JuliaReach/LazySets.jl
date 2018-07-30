# ===================================
# Approximations in the infinity norm
# ===================================

"""
    box_approximation(S::LazySet)::Hyperrectangle

Overapproximate a convex set by a tight hyperrectangle.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle.

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given set in the canonical directions, and the lengths of the sides can
be recovered from the distance among support functions in the same directions.
"""
function box_approximation(S::LazySet)::Hyperrectangle
    (c, r) = box_approximation_helper(S)
    return Hyperrectangle(c, r)
end

# special case: Hyperrectangle
box_approximation(S::Hyperrectangle) = S

# special case: other rectangle
box_approximation(S::AbstractHyperrectangle) =
    Hyperrectangle(center(S), radius_hyperrectangle(S))

"""
    interval_hull

Alias for `box_approximation`.
"""
interval_hull = box_approximation

"""
    box_approximation_symmetric(S::LazySet{N})::Hyperrectangle{N} where {N<:Real}

Overapproximate a convex set by a tight hyperrectangle centered in the origin.

### Input

- `S` -- convex set

### Output

A tight hyperrectangle centered in the origin.

### Algorithm

The center of the box is the origin, and the radius is obtained by computing the
maximum value of the support function evaluated at the canonical directions.
"""
function box_approximation_symmetric(S::LazySet{N}
                                    )::Hyperrectangle{N} where {N<:Real}
    (c, r) = box_approximation_helper(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

function box_approximation_symmetric_parallel(S::LazySet{N}
                                    )::Hyperrectangle{N} where {N<:Real}
    (c, r) = box_approximation_helper_parallel(S)
    return Hyperrectangle(zeros(N, length(c)), abs.(c) .+ r)
end

"""
    symmetric_interval_hull

Alias for `box_approximation_symmetric`.
"""
symmetric_interval_hull = box_approximation_symmetric
symmetric_interval_hull_parallel = box_approximation_symmetric_parallel

"""
    box_approximation_helper(S::LazySet)

Common code of `box_approximation` and `box_approximation_symmetric`.

### Input

- `S` -- convex set

### Output

A tuple containing the data that is needed to construct a tightly
overapproximating hyperrectangle.

- `c` -- center
- `r` -- radius

### Algorithm

The center of the hyperrectangle is obtained by averaging the support function
of the given convex set in the canonical directions.
The lengths of the sides can be recovered from the distance among support
functions in the same directions.
"""
@inline function box_approximation_helper(S::LazySet{N}) where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    n = dim(S)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    d = zeros(N, n)
    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ(d, S)
        d[i] = -one_N
        hbottom = -ρ(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        r[i] = (htop - hbottom) / 2
    end
    return c, r
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

# This function retuns the (irange,jrange) indexes assigned to this worker
function myrange(c::SharedVector{N}) where {N<:Real}
    idx = indexpids(c)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(c))
    splits = [round(Int, s) for s in linspace(0,size(c,1),nchunks+1)]
    splits[idx]+1:splits[idx+1]
end

assign_chunk!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}) where {N<:Real} = process_chunk!(c, r, S, myrange(c))

function distribute_task!(c::SharedVector{N}, r::SharedVector{N}, S::LazySet{N}) where {N<:Real}
    @sync begin
        for p in procs(c)
            @async remotecall_wait(assign_chunk!, p, c, r, S)
        end
    end
end


"""
    ballinf_approximation(S::LazySet{N})::BallInf{N} where {N<:Real}

Overapproximate a convex set by a tight ball in the infinity norm.

### Input

- `S` -- convex set

### Output

A tight ball in the infinity norm.

### Algorithm

The center and radius of the box are obtained by evaluating the support function
of the given convex set along the canonical directions.
"""
function ballinf_approximation(S::LazySet{N})::BallInf{N} where {N<:Real}
    zero_N = zero(N)
    one_N = one(N)
    n = dim(S)
    c = Vector{N}(undef, n)
    r = zero_N
    d = zeros(N, n)

    @inbounds for i in 1:n
        d[i] = one_N
        htop = ρ(d, S)
        d[i] = -one_N
        hbottom = -ρ(d, S)
        d[i] = zero_N
        c[i] = (htop + hbottom) / 2
        rcur = (htop - hbottom) / 2
        if (rcur > r)
            r = rcur
        end
    end
    return BallInf(c, r)
end
