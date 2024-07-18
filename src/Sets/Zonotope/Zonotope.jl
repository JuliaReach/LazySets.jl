"""
    Zonotope{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}} <: AbstractZonotope{N}

Type that represents a zonotope.

### Fields

- `center`     -- center of the zonotope
- `generators` -- matrix; each column is a generator of the zonotope

### Notes

Mathematically, a zonotope is defined as the set

```math
Z = \\left\\{ x ∈ ℝ^n : x = c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i ∈ [-1, 1]~~ ∀ i = 1,…, p \\right\\},
```
where ``c ∈ ℝ^n`` is its *center* and ``\\{g_i\\}_{i=1}^p``,
``g_i ∈ ℝ^n``, is the set of *generators*.
This characterization defines a zonotope as the finite Minkowski sum of line
segments.
Zonotopes can be equivalently described as the image of a unit infinity-norm
ball in ``ℝ^n`` by an affine transformation.

Zonotopes can be constructed in two different ways: either passing the
generators as a matrix, where each column represents a generator, or passing a
list of vectors, where each vector represents a generator. Below we illustrate
both ways.

### Examples

A two-dimensional zonotope with given center and matrix of generators:

```jldoctest zonotope_label
julia> Z = Zonotope([1.0, 0.0], [0.1 0.0; 0.0 0.1])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 0.0], [0.1 0.0; 0.0 0.1])

julia> dim(Z)
2

julia> center(Z)
2-element Vector{Float64}:
 1.0
 0.0

julia> genmat(Z)
2×2 Matrix{Float64}:
 0.1  0.0
 0.0  0.1
```
Here, the first vector in the `Zonotope` constructor corresponds to the center
and each column of the second argument corresponds to a generator. The functions
`center` and `genmat` respectively return the center and the generator matrix of
a zonotope.

We can collect the vertices using `vertices_list`:

```jldoctest zonotope_label
julia> vertices_list(Z)
4-element Vector{Vector{Float64}}:
 [1.1, 0.1]
 [0.9, 0.1]
 [0.9, -0.1]
 [1.1, -0.1]
```

The support vector along a given direction can be computed using `σ`
(resp. the support function can be computed using `ρ`):

```jldoctest zonotope_label
julia> σ([1.0, 1.0], Z)
2-element Vector{Float64}:
 1.1
 0.1
```

Zonotopes admit an alternative constructor that receives a list of
vectors, each vector representing a generator:

```jldoctest
julia> Z = Zonotope(ones(2), [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
Zonotope{Float64, Vector{Float64}, Matrix{Float64}}([1.0, 1.0], [1.0 0.0 1.0; 0.0 1.0 1.0])

julia> genmat(Z)
2×3 Matrix{Float64}:
 1.0  0.0  1.0
 0.0  1.0  1.0
```
"""
struct Zonotope{N,VN<:AbstractVector{N},MN<:AbstractMatrix{N}} <: AbstractZonotope{N}
    center::VN
    generators::MN

    function Zonotope(center::VN,
                      generators::MN) where {N,
                                             VN<:AbstractVector{N},
                                             MN<:AbstractMatrix{N}}
        @assert length(center) == size(generators, 1) "the dimension of the " *
                                                      "center ($(length(center))) and the generators " *
                                                      "($(size(generators, 1))) need to match"
        return new{N,VN,MN}(center, generators)
    end
end

# constructor from center and list of generators
function Zonotope(center::VN, generators_list::AbstractVector{VN}) where {VN<:AbstractVector}
    G = to_matrix(generators_list, length(center))
    return Zonotope(center, G)
end
