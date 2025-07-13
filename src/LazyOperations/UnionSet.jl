import Base: ∪

export UnionSet,
       swap

"""
    UnionSet{N, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the set union of two sets.

### Fields

- `X` -- set
- `Y` -- set

### Notes

The union of convex sets is typically not convex.
"""
struct UnionSet{N,S1<:LazySet{N},S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function UnionSet(X::LazySet{N}, Y::LazySet{N}) where {N}
        @assert dim(X) == dim(Y) "sets in a union must have the same dimension"
        return new{N,typeof(X),typeof(Y)}(X, Y)
    end
end

isoperationtype(::Type{<:UnionSet}) = true

isconvextype(::Type{<:UnionSet}) = false

function ispolyhedral(::UnionSet)
    throw(ArgumentError("this operation is not implemented"))
end

# EmptySet is the neutral element for UnionSet
@neutral(UnionSet, EmptySet)

# Universe is the absorbing element for UnionSet
@absorbing(UnionSet, Universe)

# interface for binary set operations
Base.first(U::UnionSet) = U.X
second(U::UnionSet) = U.Y
@declare_binary_operation(UnionSet)

"""
    ∪

Alias for `UnionSet`.
"""
∪(X::LazySet, Y::LazySet) = UnionSet(X, Y)

"""
    swap(cup::UnionSet)

Return a new `UnionSet` object with the arguments swapped.

### Input

- `cup` -- union of two sets

### Output

A new `UnionSet` object with the arguments swapped.
"""
function swap(cup::UnionSet)
    return UnionSet(cup.Y, cup.X)
end

concretize(cup::UnionSet) = UnionSet(concretize(first(cup)), concretize(second(cup)))

"""
    dim(cup::UnionSet)

Return the dimension of the union of two sets.

### Input

- `cup` -- union of two sets

### Output

The ambient dimension of the union of two sets.
"""
function dim(cup::UnionSet)
    return dim(cup.X)
end

"""
    σ(d::AbstractVector, cup::UnionSet; [algorithm]="support_vector")

Return a support vector of the union of two sets in a given direction.

### Input

- `d`         -- direction
- `cup`       -- union of two sets
- `algorithm` -- (optional, default: "support_vector"): the algorithm to compute
                 the support vector; if "support_vector", use the support
                 vector of each argument; if "support_function" use the support
                 function of each argument and evaluate the support vector of
                 only one of them

### Output

A support vector in the given direction.

### Algorithm

The support vector of the union of two sets ``X`` and ``Y`` can be obtained as
the vector that maximizes the support function of either ``X`` or ``Y``, i.e.,
it is sufficient to find the ``\\argmax(ρ(d, X), ρ(d, Y)])`` and evaluate its
support vector.

The default implementation, with option `algorithm="support_vector"`, computes
the support vector of ``X`` and ``Y`` and then compares the support function
using a dot product.

If the support function can be computed more efficiently, the alternative
implementation `algorithm="support_function"` can be used, which evaluates the
support function of each set directly and then calls only the support vector of
either ``X`` *or* ``Y``.
"""
function σ(d::AbstractVector, cup::UnionSet; algorithm="support_vector")
    X, Y = cup.X, cup.Y

    if algorithm == "support_vector"
        σX, σY = σ(d, X), σ(d, Y)
        return dot(d, σX) > dot(d, σY) ? σX : σY

    elseif algorithm == "support_function"
        ρX, ρY = ρ(d, X), ρ(d, Y)
        return ρX > ρY ? σ(d, X) : σ(d, Y)

    else
        error("algorithm $algorithm for the support vector of a `UnionSet` is" *
              "unknown")
    end
end

"""
    ρ(d::AbstractVector, cup::UnionSet)

Evaluate the support function of the union of two sets in a given direction.

### Input

- `d`   -- direction
- `cup` -- union of two sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function of the union of two sets ``X`` and ``Y`` evaluates to the
maximum of the support-function evaluations of ``X`` and ``Y``.
"""
@validate function ρ(d::AbstractVector, cup::UnionSet)
    X, Y = cup.X, cup.Y
    return max(ρ(d, X), ρ(d, Y))
end

"""
    an_element(cup::UnionSet)

Return some element of the union of two sets.

### Input

- `cup` -- union of two sets

### Output

An element in the union of two sets.

### Algorithm

We use `an_element` on the first non-empty wrapped set.
"""
function an_element(cup::UnionSet)
    if isempty(cup.X)
        return an_element(cup.Y)
    end
    return an_element(cup.X)
end

"""
    ∈(x::AbstractVector, cup::UnionSet)

Check whether a given point is contained in the union of two sets.

### Input

- `x`   -- point/vector
- `cup` -- union of two sets

### Output

`true` iff ``x ∈ cup``.
"""
@validate function ∈(x::AbstractVector, cup::UnionSet)
    return x ∈ cup.X || x ∈ cup.Y
end

"""
    isempty(cup::UnionSet)

Check whether the union of two sets is empty.

### Input

- `cup` -- union of two sets

### Output

`true` iff the union is empty.
"""
function isempty(cup::UnionSet)
    return isempty(cup.X) && isempty(cup.Y)
end

"""
    isbounded(cup::UnionSet)

Check whether the union of two sets is bounded.

### Input

- `cup` -- union of two sets

### Output

`true` iff the union is bounded.
"""
function isbounded(cup::UnionSet)
    return isbounded(cup.X) && isbounded(cup.Y)
end

function isboundedtype(::Type{<:UnionSet{N,S1,S2}}) where {N,S1,S2}
    return isboundedtype(S1) && isboundedtype(S2)
end

"""
    vertices_list(cup::UnionSet; [apply_convex_hull]::Bool=false,
                  [backend]=nothing)

Return the list of vertices of the union of two sets.

### Input

- `cup`               -- union of two sets
- `apply_convex_hull` -- (optional, default: `false`) if `true`, post-process
                         the vertices using a convex-hull algorithm
- `backend`           -- (optional, default: `nothing`) backend for computing
                         the convex hull (see argument `apply_convex_hull`)

### Output

The list of vertices, possibly reduced to the list of vertices of the convex
hull.
"""
function vertices_list(cup::UnionSet;
                       apply_convex_hull::Bool=false,
                       backend=nothing)
    vlist = vcat(vertices_list(cup.X), vertices_list(cup.Y))
    if apply_convex_hull
        convex_hull!(vlist; backend=backend)
    end
    return vlist
end

function linear_map(M::AbstractMatrix, cup::UnionSet)
    return UnionSet(linear_map(M, cup.X), linear_map(M, cup.Y))
end

function volume(cup::UnionSet)
    return volume(cup.X) + volume(cup.Y) - volume(Intersection(cup.X, cup.Y))
end

function convex_hull(cup::UnionSet)
    return convex_hull(cup.X, cup.Y)
end

function translate(cup::UnionSet, x::AbstractVector)
    X = translate(first(cup), x)
    Y = translate(second(cup), x)
    return UnionSet(X, Y)
end
