"""
# Extended help

    project(H::HalfSpace, block::AbstractVector{Int}; [kwargs...])

### Algorithm

If the unconstrained dimensions of `H` are a subset of the `block` variables,
the projection is applied to the normal direction of `H`.
Otherwise, the projection results in the universal set.

The latter can be seen as follows.
Without loss of generality consider projecting out a single and constrained
dimension ``xₖ`` (projecting out multiple dimensions can be modeled by
repeatedly projecting out one dimension).
We can write the projection as an existentially quantified linear constraint:

```math
    ∃xₖ: a₁x₁ + … + aₖxₖ + … + aₙxₙ ≤ b
```

Since ``aₖ ≠ 0``, there is always a value for ``xₖ`` that satisfies the
constraint for any valuation of the other variables.

### Examples

Consider the half-space ``x + y + 0⋅z ≤ 1``, whose ambient dimension is `3`.
The (trivial) projection in the three dimensions using the block of variables
`[1, 2, 3]` is:

```jldoctest project_halfspace
julia> H = HalfSpace([1.0, 1.0, 0.0], 1.0)
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, 0.0], 1.0)

julia> project(H, [1, 2, 3])
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0, 0.0], 1.0)
```

Projecting along dimensions `1` and `2` only:

```jldoctest project_halfspace
julia> project(H, [1, 2])
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```

For convenience, one can use `project(H, constrained_dimensions(H))` to return
the half-space projected on the dimensions where it is constrained:

```jldoctest project_halfspace
julia> project(H, constrained_dimensions(H))
HalfSpace{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```

If a constrained dimension is projected, we get the universal set of the
dimension corresponding to the projection.

```jldoctest project_halfspace
julia> project(H, [1, 3])
Universe{Float64}(2)

julia> project(H, [1])
Universe{Float64}(1)
```
"""
@validate function project(H::HalfSpace, block::AbstractVector{Int}; kwargs...)
    if constrained_dimensions(H) ⊆ block
        return HalfSpace(H.a[block], H.b)
    else
        N = eltype(H)
        return Universe{N}(length(block))
    end
end
