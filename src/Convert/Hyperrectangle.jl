"""
    convert(::Type{Hyperrectangle}, cpa::CartesianProductArray{N, HN})
        where {N, HN<:AbstractHyperrectangle}

Convert the Cartesian product of a finite number of hyperrectangular sets to
a single hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product array of hyperrectangular set

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N,HN}) where {N,HN<:AbstractHyperrectangle}
    n = dim(cpa)
    c = Vector{N}(undef, n)
    r = Vector{N}(undef, n)
    i = 1
    @inbounds for block_set in cpa
        j = i + dim(block_set) - 1
        c[i:j] = center(block_set)
        r[i:j] = radius_hyperrectangle(block_set)
        i = j + 1
    end
    return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle}, cp::CartesianProduct{N, HN1, HN2})
        where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}

Convert the Cartesian product of two hyperrectangular sets to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product of two hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

The result is obtained by concatenating the center and radius of each
hyperrectangle. This implementation uses the `center` and
`radius_hyperrectangle` methods.
"""
function convert(::Type{Hyperrectangle},
                 cp::CartesianProduct{N,HN1,HN2}) where {N,HN1<:AbstractHyperrectangle,
                                                         HN2<:AbstractHyperrectangle}
    X, Y = first(cp), second(cp)
    c = vcat(center(X), center(Y))
    r = vcat(radius_hyperrectangle(X), radius_hyperrectangle(Y))
    return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, IN}) where {N, IN<:Interval}

Convert the Cartesian product of a finite number of intervals to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `S`              -- Cartesian product array of intervals

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `min` and `max` methods of `Interval` to reduce
the allocations and improve performance (see LazySets#1143).
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N,IN}) where {N,IN<:Interval}
    # since the sets are intervals, the dimension of cpa is its length
    n = length(array(cpa))
    l = Vector{N}(undef, n)
    h = Vector{N}(undef, n)
    @inbounds for (i, Ii) in enumerate(array(cpa))
        l[i] = min(Ii)
        h[i] = max(Ii)
    end
    return Hyperrectangle(; low=l, high=h)
end

"""
    convert(::Type{Hyperrectangle}, r::Rectification{N, AH})
        where {N, AH<:AbstractHyperrectangle}

Convert a rectification of a hyperrectangle to a hyperrectangle.

### Input

- `Hyperrectangle` -- target type
- `r`              -- rectification of a hyperrectangle

### Output

A `Hyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 r::Rectification{N,AH}) where {N,AH<:AbstractHyperrectangle}
    return rectify(r.X)
end

function convert(::Type{Hyperrectangle}, Z::AbstractZonotope{N}) where {N}
    c = center(Z)
    n = length(c)
    r = zeros(N, n)
    @inbounds for cj in generators(Z)
        i = findfirst(!=(zero(N)), cj)
        if isnothing(i)
            continue
        end
        @assert isnothing(findfirst(!=(zero(N)), @view cj[(i + 1):end])) "the zonotope " *
                                                                         "is not hyperrectangular"
        r[i] += cj[i]  # `+` because to allow for multiple generators in dimension i
    end
    return Hyperrectangle(c, r)
end
