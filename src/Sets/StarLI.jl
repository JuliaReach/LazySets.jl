export StarLI,
       center,
       basis,
       predicate,
       intersection!,
       dim

"""
    StarLI{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, S1<:LazySet{N}, S2<:LazySet{N}, IT<:Intersection{N, S1, S2}} <: AbstractPolyhedron{N}

Generalized star set with a lazy predicate, i.e.

```math
X = \\{x ∈ \\mathbb{R}^n : x = x₀ + \\sum_{i=1}^m α_i v_i,~~\\textrm{s.t. } P(α) = ⊤ \\},
```
where ``x₀ ∈ \\mathbb{R}^n`` is the center, the ``m`` vectors ``v₁, …, vₘ`` form
the basis of the star set, and the combination factors
``α = (α₁, …, αₘ) ∈ \\mathbb{R}^m`` are the predicates' decision variables,
i.e. ``P : α ∈ \\mathbb{R}^m → \\{⊤, ⊥\\}`` where the lazy predicate is such that
``P(α) = ⊤`` if and only if ``α ∈ Y ∩ Z`` for given sets ``Y`` and ``Z``.

### Fields

- `c` -- vector that represents the center
- `V` -- matrix where each column corresponds to a basis vector
- `P` -- lazy intersection that represents the predicate
"""
struct StarLI{N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, S1<:LazySet{N}, S2<:LazySet{N}, IT<:Intersection{N, S1, S2}} <: AbstractPolyhedron{N}
    c::VN # center
    V::MN # basis
    P::IT # predicate

    # default constructor with size checks
    function StarLI(c::VN, V::MN, P::IT) where {N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, S1<:LazySet{N}, S2<:LazySet{N}, IT<:Intersection{N, S1, S2}}
        @assert length(c) == size(V, 1) "the center of the basis vectors should be compatible, " *
                                        "but they are of length $(length(c)) and $(size(V, 1)) respectively"

        @assert dim(P) == size(V, 2) "the number of basis vectors should be compatible " *
                                     "with the predicates' dimension, but they are $(size(V, 2)) and $(dim(P)) respectively"

        return new{N, VN, MN, S1, S2, IT}(c, V, P)
    end
end

# constructor from center and list of generators
function StarLI(c::VN, vlist::AbstractVector{VN}, P::IT) where {N, VN<:AbstractVector{N}, IT<:Intersection{N}}
    V = to_matrix(vlist, length(c))
    return Star(c, V, P)
end

# constructor from a zonotope (no additional constraints)
# default constructor with size checks
function StarLI(c::VN, V::MN, Z::ZT) where {N, VN<:AbstractVector{N}, MN<:AbstractMatrix{N}, ZT<:AbstractZonotope{N}}
    Pempty = HPolyhedron{N, Vector{N}}()
    return StarLI(c, V, Intersection(Z, Pempty, check=false))
end

isoperationtype(::Type{<:StarLI}) = false
isconvextype(::Type{<:StarLI}) = true

# =================================================
# Star set with lazy constraints getter functions
# =================================================

center(X::StarLI) = X.c

basis(X::StarLI) = X.V

predicate(X::StarLI) = X.P

# =============================
# LazySets interface functions
# =============================

function dim(X::StarLI)
    return length(X.c)
end

function ρ(d::AbstractVector, X::StarLI)
    P = concretize(predicate(X))
    μ = ρ(_At_mul_B(basis(X), d), P)
    return μ + dot(d, center(X))
end

function ρ_upper_bound(d::AbstractVector, X::StarLI)
    P = predicate(X)
    a = ρ(_At_mul_B(basis(X), d), P.X)
    b = ρ(_At_mul_B(basis(X), d), P.Y)
    μ = min(a, b)
    return μ + dot(d, center(X))
end

# TODO
#function box_approximation(X::StarLI)
    # use ρ_upper_bound
#end

function σ(d::AbstractVector, X::StarLI)
    # exact
    A = basis(X)
    P = concretize(predicate(X))
    return A * σ(_At_mul_B(A, d), P) + center(X)
end

function an_element(X::StarLI)
    return basis(X) * an_element(predicate(X)) + center(X)
end

function isempty(X::StarLI)
    return isempty(predicate(X))
end

function ∈(x::AbstractVector, X::StarLI)
    return basis(X) \ (x - center(X)) ∈ predicate(X)
end

function vertices_list(X::StarLI; apply_convex_hull::Bool=true)
    P = predicate(X)
    Pint = intersection(P.X, P.Y)
    S = Star(X.c, X.V, Pint)
    am = convert(STAR, S)
    return vertices_list(am, apply_convex_hull=apply_convex_hull)
end

function constraints_list(X::StarLI)
    P = predicate(X)
    Pint = intersection(P.X, P.Y)
    S = Star(X.c, X.V, Pint)
    am = convert(STAR, X)
    return constraints_list(am)
end

function linear_map(M::AbstractMatrix, X::StarLI)
    c′ = M * X.c
    V′ = M * X.V
    P′ = X.P
    return StarLI(c′, V′, P′)
end

function affine_map(M::AbstractMatrix, X::StarLI, v::AbstractVector)
    c′ = M * X.c + v
    V′ = M * X.V
    P′ = X.P
    return Star(c′, V′, P′)
end

function rand(::Type{StarLI};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing)
    rng = reseed(rng, seed)
    Z = rand(Zonotope, N=N, dim=dim, rng=rng, seed=seed)
    X = convert(StarLI, Z)
    return X
end

# ============================
# Concrete intersection
# ============================

function intersection!(X::StarLI{N, VN, MN, S1, S2}, H::HalfSpace) where {N, VN, MN, S1<:AbstractZonotope{N}, S2<:HPolyhedron{N}}
    P = predicate(X)
    _intersection_star!(center(X), basis(X), P.Y, H)
    return X
end

function intersection(X::StarLI{N, VN, MN, S1, S2}, H::HalfSpace) where {N, VN, MN, S1<:AbstractZonotope{N}, S2<:HPolyhedron{N}}
    c = center(X)
    V = basis(X)
    P = predicate(X)

    Pnew = copy(P)
    Xnew = StarLI(c, V, Pnew)
    return intersection!(Xnew, H)
end

# symmetric methods
intersection(H::HalfSpace, X::StarLI) = intersection(X, H)
