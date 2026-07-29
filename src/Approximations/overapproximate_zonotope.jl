"""
    overapproximate(Z::AbstractZonotope, ::Type{<:Zonotope}, r::Real)

Reduce the order of a zonotopic set.

### Input

- `Z`        -- zonotopic set
- `Zonotope` -- target set type
- `r`        -- desired order

### Output

A new zonotope with `r` generators, if possible.

### Algorithm

This method falls back to `reduce_order` with the default algorithm.
"""
function overapproximate(Z::AbstractZonotope, ::Type{<:Zonotope}, r::Real)
    return reduce_order(Z, r)
end

"""
    overapproximate(X::ConvexHull{N, <:AbstractZonotope, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the convex hull of two zonotopic sets.

### Input

- `X`         -- convex hull of two zonotopic sets
- `Zonotope`  -- target set type
- `algorithm` -- (optional; default: `"mean"`) choice of algorithm; possible
                 values are `"mean"` and `"join"`

### Output

A zonotope ``Z`` such that ``X ⊆ Z``.

### Algorithm

The algorithm can be controlled by the parameter `algorithm`.
Note that the results of the two implemented algorithms are generally
incomparable.

##### 'mean' method

If `algorithm == "mean"`, we choose the method proposed in [Girard05](@citet).
The convex hull of two zonotopic sets ``Z₁`` and ``Z₂`` of the same order,
which we write

```math
Z_j = ⟨c^{(j)}, g^{(j)}_1, …, g^{(j)}_p⟩
```
for ``j = 1, 2``, can be overapproximated as follows:

```math
CH(Z_1, Z_2) ⊆ \\frac{1}{2}⟨c^{(1)}+c^{(2)}, g^{(1)}_1+g^{(2)}_1, …,
g^{(1)}_p+g^{(2)}_p, c^{(1)}-c^{(2)}, g^{(1)}_1-g^{(2)}_1, …, g^{(1)}_p-g^{(2)}_p⟩.
```

If the zonotope order is not the same, this algorithm calls
`reduce_order` to reduce the order to the minimum of the arguments.

It should be noted that the output zonotope is not necessarily the minimal
enclosing zonotope, which is in general expensive to obtain in high dimensions.
This is further investigated in [GuibasNZ03](@citet).

##### 'join' method

If `algorithm == "join"`, we choose the method proposed in [GhorbalGP09; Definition 1](@citet).
The convex hull ``X`` of two zonotopic sets ``Z₁`` and ``Z₂`` is
overapproximated by a zonotope ``Z₃`` such that the box approximation of ``X``
is identical with the box approximation of ``Z₃``.
Let ``□(X)`` denote the box approximation of ``X``.
The center of ``Z₃`` is the center of ``□(X)``.

The generator construction consists of two phases.
In the first phase, we construct generators ``g`` as a combination of one
generator from ``Z₁``, say, ``g₁``, with another generator from ``Z₂``, say,
``g₂``.
The entry of ``g`` in the ``i``-th dimension is given as

```math
    g[i] = \\argmin_{\\min(g₁[i], g₂[i]) ≤ x ≤ \\max(g₁[i], g₂[i])} |x|.
```

If ``g`` is the zero vector, it can be omitted.

In the second phase, we construct another generator for each dimension.
These generators are scaled unit vectors.
The following formula defines the sum of all those generators.

```math
    \\sup(□(X)) - c - ∑_g |g|
```

where ``c`` is the center of the new zonotope and the ``g``s are the generators
constructed in the first phase.
"""
function overapproximate(X::ConvexHull{N,<:AbstractZonotope,<:AbstractZonotope},
                         ::Type{<:Zonotope}; algorithm="mean") where {N}
    return _overapproximate_union_zonotope(X, algorithm)
end

function overapproximate(X::UnionSet{N,<:AbstractZonotope,<:AbstractZonotope},
                         ::Type{<:Zonotope}; algorithm="mean") where {N}
    return _overapproximate_union_zonotope(X, algorithm)
end

function _overapproximate_union_zonotope(X::LazySet, algorithm)
    # execute specific algorithm
    if algorithm == "mean"
        return _overapproximate_union_zonotope_G05(X)
    elseif algorithm == "join"
        return _overapproximate_union_zonotope_GGP09(X)
    else
        throw(ArgumentError("algorithm $algorithm unknown"))
    end
end

function _overapproximate_union_zonotope_G05(X::LazySet{N}) where {N}
    # reduce to the same order if possible
    Z1, Z2 = first(X), second(X)
    m1, m2 = ngens(Z1), ngens(Z2)
    if m1 < m2
        Z2 = reduce_order(Z2, max(1, m1))
    elseif m1 > m2
        Z1 = reduce_order(Z1, max(1, m2))
    end

    if order(Z2) > order(Z1)
        Z1, Z2 = Z2, Z1
    end

    c1 = center(Z1)
    c2 = center(Z2)
    G1 = genmat(Z1)
    G2 = genmat(Z2)
    c = (c1 + c2) / N(2)

    # the case of equal order is treated separately to avoid a slicing
    # (this creates a copy)
    if order(Z1) == order(Z2)
        G = hcat(G1 .+ G2,
                 c1 - c2,
                 G1 .- G2) / N(2)
    else
        G = hcat((G1[:, 1:ngens(Z2)] .+ G2) / N(2),
                 (c1 - c2) / N(2),
                 (G1[:, 1:ngens(Z2)] .- G2) / N(2),
                 G1[:, (ngens(Z2) + 1):end])
    end
    Z = Zonotope(c, G)
    return remove_zero_generators(Z)
end

function _overapproximate_union_zonotope_GGP09(X::LazySet{N}) where {N}
    Z1, Z2 = first(X), second(X)
    m = min(ngens(Z1), ngens(Z2))
    G1, G2 = genmat(Z1), genmat(Z2)
    n = dim(Z1)
    box = box_approximation(X)

    # new center: mid point of box approximation
    c = center(box)

    # the k-th new generator is a simple combination of the old k-th generators
    G = Vector{Vector{N}}()
    sizehint!(G, m + n)
    g_sum = zeros(N, n)
    @inbounds for j in 1:m
        g1, g2 = G1[:, j], G2[:, j]
        g = Vector{N}(undef, n)
        for i in 1:n
            gi_min, gi_max = g1[i] < g2[i] ? (g1[i], g2[i]) : (g2[i], g1[i])
            if gi_min <= zero(N) && gi_max >= zero(N)
                g[i] = zero(N)
            elseif abs(gi_min) <= abs(gi_max)
                g[i] = gi_min
            else
                g[i] = gi_max
            end
        end
        if !iszero(g)
            push!(G, g)
            g_sum += abs.(g)
        end
    end

    # one more new generator (a scaled unit vector) for every dimension
    g_total = high(box) - c - g_sum
    @inbounds for i in 1:n
        if !iszero(g_total[i])
            g = SingleEntryVector(i, n, g_total[i])
            push!(G, g)
        end
    end

    return Zonotope(c, G)
end

# Given center, (dependent) generator matrix and exponent matrix of a (simple)
# sparse polynomial zonotope, compute the new center and generator matrix of
# its zonotope overapproximation. This method assumes that the parameter domain
# is [-1, 1]ᵖ.
function _zonotope_overapprox!(c, G, E)
    @inbounds for (j, g) in enumerate(eachcol(G))
        if all(iseven, E[:, j])
            c .+= g / 2
            G[:, j] ./= 2
        end
    end
    return c, G
end

function _zonotope_overapprox(c, G, E)
    return _zonotope_overapprox!(copy(c), copy(G), E)
end

"""
    overapproximate(P::AbstractSparsePolynomialZonotope, ::Type{<:Zonotope})

Overapproximate a sparse polynomial zonotope with a zonotope.

### Input

- `P`        -- sparse polynomial zonotope
- `Zonotope` -- target set type

### Output

A zonotope.

### Algorithm

This method implements [Kochdumper21a; Proposition 3.1.14](@citet).
"""
function overapproximate(P::AbstractSparsePolynomialZonotope, ::Type{<:Zonotope})
    cnew, Gnew = _zonotope_overapprox(center(P), genmat_dep(P), expmat(P))
    if ngens_indep(P) > 0
        Gnew = remove_redundant_generators(hcat(Gnew, genmat_indep(P)))
    end
    return Zonotope(cnew, Gnew)
end

"""
    overapproximate(P::AbstractSparsePolynomialZonotope, ::Type{<:Zonotope},
                    dom::AbstractVector{<:IA.Interval})

Overapproximate a sparse polynomial zonotope over the parameter domain `dom`
with a zonotope.

### Input

- `P`         -- sparse polynomial zonotope
- `Zonotope`  -- target set type
- `dom`       -- parameter domain, which should be a subset of `[-1, 1]^q`,
                 where `q = nparams(P)`

### Output

A zonotope.
"""
function overapproximate(P::AbstractSparsePolynomialZonotope, ::Type{<:Zonotope},
                         dom::AbstractVector{<:IA.Interval})
    @assert length(dom) == nparams(P) &&
            all(IA.issubset_interval(vi, IA.interval(-1, 1)) for vi in dom) "dom should be a " *
                                                                            "subset of [-1, 1]^q"

    # handle dependent generators
    G = genmat_dep(P)
    E = expmat(P)
    cnew = copy(center(P))
    Gnew = similar(G)
    @inbounds for (j, g) in enumerate(eachcol(G))
        # monomial value over the domain
        α = IA.interval(1, 1)
        for (i, vi) in enumerate(dom)
            α *= vi^E[i, j]
        end
        m, r = IA.midradius(α)
        cnew .+= m * g
        Gnew[:, j] .= r * g
    end

    # handle independent generators
    if ngens_indep(P) > 0
        Gnew = remove_redundant_generators(hcat(Gnew, genmat_indep(P)))
    end

    return Zonotope(cnew, Gnew)
end

"""
    overapproximate(P::DensePolynomialZonotope, ::Type{<:Zonotope})

Overapproximate a polynomial zonotope with a zonotope.

### Input

- `P`        -- polynomial zonotope
- `Zonotope` -- target set type

### Output

A zonotope.

### Algorithm

This method implements [Althoff13; Proposition 1](@citet).
"""
function overapproximate(P::DensePolynomialZonotope, ::Type{<:Zonotope})
    η = polynomial_order(P)
    cnew = center(P) + 1 / 2 * vec(sum(i -> sum(P.E[2i]; dims=2), 1:floor(Int, η / 2)))
    Gnew = hcat([iseven(i) ? 1 / 2 * P.E[i] : P.E[i] for i in 1:η]..., P.F..., P.G)
    return Zonotope(cnew, Gnew)
end

"""
    overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                    dir::Union{AbstractDirections, Type{<:AbstractDirections}};
                    [algorithm]="vrep", kwargs...)

Overapproximate a set with a zonotope.

### Input

- `X`         -- set
- `Zonotope`  -- target set type
- `dir`       -- directions used for the generators
- `algorithm` -- (optional, default: `"vrep"`) algorithm used to compute the
                 overapproximation
- `kwargs`    -- further algorithm-specific options

### Output

A zonotope that overapproximates `X` and uses at most the generator directions
provided in `dir` (redundant directions will be ignored).

### Notes

Two algorithms are available:

- `"vrep"` -- Overapproximate a polytopic set with a zonotope of minimal total
  generator sum using only generators in the given directions. Under this
  constraint, the zonotope has the minimal sum of generator vectors. See the
  docstring of [`_overapproximate_zonotope_vrep`](@ref) for further details.

- `"cpa"` -- Overapproximate a polytopic set with a zonotope using a Cartesian
  decomposition into two-dimensional blocks. See the docstring of
  [`_overapproximate_zonotope_cpa`](@ref) for further details.
"""
function overapproximate(X::LazySet, ZT::Type{<:Zonotope},
                         dir::Union{AbstractDirections,Type{<:AbstractDirections}};
                         algorithm="vrep", kwargs...)
    if algorithm == "vrep"
        return _overapproximate_zonotope_vrep(X, dir, kwargs...)
    elseif algorithm == "cpa"
        cpa = _overapproximate_zonotope_cpa(X, dir)
        return convert(Zonotope, cpa)
    else
        throw(ArgumentError("algorithm $algorithm unknown"))
    end
end

# disambiguation
function overapproximate(∅::EmptySet, ::Type{<:Zonotope},
                         dir::AbstractDirections;
                         algorithm="vrep", kwargs...)
    return ∅
end
function overapproximate(∅::EmptySet, ::Type{<:Zonotope},
                         dir::Type{<:AbstractDirections};
                         algorithm="vrep", kwargs...)
    return ∅
end

"""
    _overapproximate_zonotope_vrep(X::LazySet{N},
                                   dir::AbstractDirections;
                                   solver=default_lp_solver(N)) where {N}

Overapproximate a polytopic set with a zonotope of minimal total generator sum
using only generators in the given directions.

### Input

- `X`      -- polytopic set
- `dir`    -- directions used for the generators
- `solver` -- (optional, default: `default_lp_solver(N)`) the backend used to
              solve the linear program

### Output

A zonotope that overapproximates `X` and uses at most the directions provided in
`dir` (redundant directions will be ignored).
Under this constraint, the zonotope has the minimal sum of generator vectors.

### Notes

The algorithm only requires one representative of each generator direction and
their additive inverse (e.g., only one of `[1, 0]` and `[-1, 0]`) and assumes
that the directions are normalized.
We preprocess the directions in that respect.

### Algorithm

We solve a linear program parametric in the vertices ``v_j`` of `X` and the
directions ``d_k`` in `dir` presented in Section 4.2 in [GuibasNZ03](@citet),
adapting the notation to the one used in this library.

```math
    \\min ∑_{k=1}^l α_k \\
    s.t. \\
    c + ∑_{k=1}^l b_{kj} * d_k = v_j \\quad ∀ j \\
    -α_k ≤ b_{kj} ≤ α_k \\quad ∀ k, j \\
    α_k ≥ 0 \\quad ∀ k
```

The resulting zonotope has center `c` and generators `α_k · d_k`.

Note that the first type of side constraints is vector-based and that the
nonnegativity constraints (last type) are not stated explicitly in [GuibasNZ03](@cite).
"""
function _overapproximate_zonotope_vrep(X::LazySet{N},
                                        dir::AbstractDirections;
                                        solver=default_lp_solver(N)) where {N}
    @assert ispolytopic(X) "this algorithm requires a polytope"

    # TODO "normalization" here involves two steps: removing opposite directions
    # and normalizing the direction vector
    # for the latter we can use the normalization information from dispatch on
    # the directions (some typical such directions are normalized by definition)
    function normalize_directions(dir::AD) where {AD}
        dirs = Vector{Vector{N}}()
        for d in dir
            d_normalized = normalize(d)
            if (d_normalized ∉ dirs) && (-d_normalized ∉ dirs)
                push!(dirs, d_normalized)
            end
        end
        return dirs
    end

    n = dim(X)
    dirs = normalize_directions(dir)
    l = length(dirs)
    V = vertices_list(X)
    m = length(V)
    nvariables = l + n + l * m  # l 'α_k' + n 'p[i]' + lm 'b_lj'
    nconstraints = n * m + 2 * l * m  # nm 'p_j' + 2lm 'b_kj'

    obj = vcat(ones(l), zeros(nvariables - l))
    A = spzeros(N, nconstraints, nvariables)
    sense = vcat(fill('=', n * m), fill('<', 2 * l * m))
    b = zeros(N, nconstraints)
    r = 1
    # constraints p[i] + ∑_k v_j*b_kj = p_j
    col_offset_p = l
    col_offset_b = l + n + 1
    for (j, vj) in enumerate(V)
        for i in 1:n
            A[r, col_offset_p + i] = N(1)  # p[i]
            for (k, dk) in enumerate(dirs)
                A[r, col_offset_b + (k - 1) * m] = dk[i]
            end
            b[r] = vj[i]
            r += 1
        end
        col_offset_b += 1
    end
    # constraints -α_k ± b_kj <= 0
    col_offset_b = l + n
    for k in 1:l
        for j in 1:m
            for sign_b in N[1, -1]
                A[r, k] = N(-1)  # - α_k
                A[r, col_offset_b + j] = sign_b  # ± b_kj
                r += 1
            end
        end
        col_offset_b += m
    end
    lbounds = vcat(zeros(N, l), fill(N(-Inf), nvariables - l))
    ubounds = fill(Inf, nvariables)
    lp = linprog(obj, A, sense, b, lbounds, ubounds, solver)

    if !is_lp_optimal(lp.status)
        throw(ArgumentError("got unexpected status from LP solver: $(lp.status)"))
    end
    c = lp.sol[(col_offset_p + 1):(col_offset_p + n)]
    ck = lp.sol[1:l]
    G = Matrix{N}(undef, n, l)
    for (j, dj) in enumerate(dirs)
        G[:, j] = ck[j] * dj
    end
    return Zonotope(c, G)
end

# overload on direction type
function _overapproximate_zonotope_vrep(X::LazySet{N},
                                        dir::Type{<:AbstractDirections};
                                        solver=default_lp_solver(N)) where {N}
    return _overapproximate_zonotope_vrep(X, _get_directions(dir, dim(X)); solver=solver)
end

"""
    _overapproximate_zonotope_cpa(X::LazySet, dir::Type{<:AbstractDirections})

Overapproximate a polytopic set with a zonotope using Cartesian decomposition.

### Input

- `X`   -- polytopic set
- `dir` -- directions used for the generators

### Output

A zonotope that overapproximates `X`.

### Notes

The algorithm decomposes `X` into 2D sets and overapproximates those sets with
zonotopes, and finally converts the Cartesian product of the sets to a zonotope.

### Algorithm

The implementation is based on [LeGuernic09; Section 8.2.4](@citet).
"""
function _overapproximate_zonotope_cpa(X::LazySet, dir::AbstractDirections)
    n = dim(X)
    if dim(dir) != 2
        # try to convert to 2D directions
        dir = _get_directions(typeof(dir), 2)
    end
    # overapproximate 2D blocks
    if n > 1
        πX_2D = [project(X, [i, i + 1]) for i in 1:2:(n - 1)]
        Z_2D = [_overapproximate_zonotope_vrep(poly, dir) for poly in πX_2D]
    end

    if iseven(n)
        out = Z_2D
    else
        # odd case projects onto an interval
        πX_n = project(X, [n])
        πX_n = overapproximate(πX_n, Interval)
        πX_n = convert(Zonotope, πX_n)

        if n == 1
            out = [πX_n]
        else
            out = vcat(Z_2D, πX_n)
        end
    end
    return CartesianProductArray(out)
end

# overload on direction type
function _overapproximate_zonotope_cpa(X::LazySet, dir::Type{<:AbstractDirections})
    return _overapproximate_zonotope_cpa(X, _get_directions(dir, 2))
end

"""
    overapproximate(r::Rectification{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the rectification of a zonotopic set with a zonotope.

### Input

- `r`        -- lazy rectification of a zonotopic set
- `Zonotope` -- target set type

### Output

A zonotope overapproximation of the set obtained by rectifying `Z`.

### Algorithm

This method implements [SinghGMPV18; Theorem 3.1](@citet).
"""
function overapproximate(r::Rectification{N,<:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    Z = set(r)
    c = copy(center(Z))
    G = copy(genmat(Z))
    n, m = size(G)
    row_idx = Vector{Int}()
    μ_idx = Vector{N}()

    @inbounds for i in 1:n
        lx, ux = low(Z, i), high(Z, i)
        if !_leq(lx, zero(N))
            nothing
        elseif _leq(ux, zero(N))
            c[i] = zero(N)
            for j in 1:m
                G[i, j] = zero(N)
            end
        else
            λ = ux / (ux - lx)
            μ = -λ * lx / 2
            c[i] = c[i] * λ + μ
            for j in 1:m
                G[i, j] = G[i, j] * λ
            end
            push!(row_idx, i)
            push!(μ_idx, μ)
        end
    end

    q = length(row_idx)
    if q >= 1
        Gnew = zeros(N, n, q)
        j = 1
        @inbounds for i in row_idx
            Gnew[i, j] = μ_idx[j]
            j += 1
        end
        Gout = hcat(G, Gnew)
    else
        Gout = G
    end

    return Zonotope(c, remove_zero_columns(Gout))
end

"""
    overapproximate(CHA::ConvexHullArray{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate the convex hull of a list of zonotopic sets with a zonotope.

### Input

- `CHA`      -- convex hull array of zonotopic sets
- `Zonotope` -- target set type

### Output

A zonotope overapproximation of the convex hull array of zonotopic sets.

### Algorithm

This method iteratively applies the overapproximation algorithm to the
convex hull of two zonotopic sets from the given array of zonotopic sets.
"""
function overapproximate(X::ConvexHullArray{N,<:AbstractZonotope}, ::Type{<:Zonotope}) where {N}
    return _overapproximate_union_zonotope(array(X))
end

function overapproximate(X::UnionSetArray{N,<:AbstractZonotope}, ::Type{<:Zonotope}) where {N}
    return _overapproximate_union_zonotope(array(X))
end

function _overapproximate_union_zonotope(arr::AbstractVector)
    n = length(arr)
    @assert n > 0 "cannot overapproximate an empty array set"
    @inbounds if n == 1
        return convert(Zonotope, arr[1])
    else
        Zaux = overapproximate(UnionSet(arr[1], arr[2]), Zonotope)
        for k in 3:n
            Zaux = overapproximate(UnionSet(Zaux, arr[k]), Zonotope)
        end
        return Zaux
    end
end

"""
    overapproximate(QM::QuadraticMap{N, <:AbstractZonotope},
                    ::Type{<:Zonotope}) where {N}

Overapproximate a quadratic map of a zonotopic set with a zonotope.

### Input

- `QM`       -- quadratic map of a zonotopic set
- `Zonotope` -- target set type

### Output

A zonotope overapproximating the quadratic map of a zonotopic set.

### Notes

Mathematically, a quadratic map of a zonotope with matrices ``Q`` is defined as:

```math
    Z_Q = \\right\\{ λ \\mid λ_i = x^T Q\\^{(i)} x,~i = 1, …, n,~x ∈ Z \\left\\}
```

### Algorithm

This method implements [AlthoffK12; Lemma 1](@citet).
"""
function overapproximate(QM::QuadraticMap{N,<:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    Z = QM.X
    G = genmat(Z)
    c = center(Z)
    n, p = size(G)
    h = Matrix{N}(undef, n, binomial(p + 2, 2) - 1)
    d = Vector{N}(undef, n)
    g(x) = view(G, :, x)
    cᵀ = c'
    for (i, Qᵢ) in enumerate(QM.Q)
        cᵀQᵢ = cᵀ * Qᵢ
        Qᵢc = Qᵢ * c
        aux = zero(N)
        for j in 1:p
            aux += g(j)' * Qᵢ * g(j)
            h[i, j] = cᵀQᵢ * g(j) + g(j)' * Qᵢc
            h[i, p + j] = g(j)' * Qᵢ * g(j) / 2
        end
        d[i] = cᵀQᵢ * c + aux / 2
        l = 0
        for j in 1:(p - 1)
            gjᵀQᵢ = g(j)' * Qᵢ
            Qᵢgj = Qᵢ * g(j)
            for k in (j + 1):p
                l += 1
                h[i, 2p + l] = gjᵀQᵢ * g(k) + g(k)' * Qᵢgj
            end
        end
    end
    return Zonotope(d, remove_zero_columns(h))
end

"""
    overapproximate(X::Intersection{N, <:AbstractZonotope, <:Hyperplane},
                    ::Type{<:Zonotope})

Overapproximate the intersection of a zonotopic set and a hyperplane with a
zonotope.

### Input

- `X`        -- intersection of a zonotopic set and a hyperplane
- `Zonotope` -- target set type

### Output

A zonotope overapproximating the intersection.

### Algorithm

This method implements [MaigaRTC14; Algorithm 3](@citet).
"""
function overapproximate(X::Intersection{N,<:AbstractZonotope,<:Hyperplane},
                         ::Type{<:Zonotope}) where {N}
    return _overapproximate_zonotope_hyperplane(first(X), second(X))
end

# symmetric method
function overapproximate(X::Intersection{N,<:Hyperplane,<:AbstractZonotope},
                         ::Type{<:Zonotope}) where {N}
    return _overapproximate_zonotope_hyperplane(second(X), first(X))
end

function _overapproximate_zonotope_hyperplane(Z::AbstractZonotope, H::Hyperplane)
    c, G = center(Z), genmat(Z)
    a, b = H.a, H.b

    s = G' * a
    d = b - dot(a, c)
    V0 = nullspace(s')

    cs = s * d / dot(s, s)
    Gs = V0 * V0'

    c_cap = c + G * cs
    G_cap = G * Gs
    return Zonotope(c_cap, G_cap)
end

function overapproximate(X::Intersection{N,<:AbstractZonotope,<:HalfSpace{N,<:SingleEntryVector}},
                         ::Type{Zonotope}) where {N}
    return _overapproximate_zonotope_halfspace_ICP(first(X), second(X))
end

# symmetric method
function overapproximate(X::Intersection{N,<:HalfSpace{N,<:SingleEntryVector},<:AbstractZonotope},
                         ::Type{Zonotope}) where {N}
    return _overapproximate_zonotope_halfspace_ICP(second(X), first(X))
end

# method for a single constraint of the form x_i <= b based on ICP
# - zonotope = G * []^p + c, where []^p is the p-dimensional unit cube
# - let G_i be the i-th row of G
# - build an ICP contractor for the constraint G_i^T x <= b - c
# - contract the domain of []^p under this constraint to D
# - the resulting zonotope is G * D + c
function _overapproximate_zonotope_halfspace_ICP(Z::AbstractZonotope{N},
                                                 H::HalfSpace{N,SingleEntryVector{N}}) where {N}
    c = center(Z)
    G = genmat(Z)
    p = size(G, 2)
    d = H.a.i

    a = G[d, :]
    io = IOBuffer()
    v = H.a.v
    negate = v < zero(N)
    if negate
        v = -v
    end
    vars = ""
    first = true
    for (i, ai) in enumerate(a)
        if first
            first = false
        else
            vars *= ", "
            if ai >= 0
                write(io, "+")
            end
        end
        write(io, "$(v * ai)*x$(i)")
        vars *= "x$(i)"
    end
    if negate
        write(io, ">=$(-H.b + c[d])")
    else
        write(io, "<=$(H.b - c[d])")
    end
    e = Meta.parse(String(take!(io)))  # NOTE: this is an internal function
    X = fill(IA.interval(-1, 1), p)
    newD = _contract_zonotope_halfspace_ICP(e, X, vars)
    if isempty(newD)
        return EmptySet{N}(length(c))
    end
    return affine_map(G, convert(Hyperrectangle, newD), c)
end

# see ext/IntervalConstraintProgrammingExt.jl
function _contract_zonotope_halfspace_ICP(e, X, vars_string)
    mod = Base.get_extension(@__MODULE__, :IntervalConstraintProgrammingExt)
    require(mod, :IntervalConstraintProgramming; fun_name="overapproximate")
    error()
end

"""
	overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                  MAT<:MatrixZonotope{NM}}

Overapproximate the linear map of a zonotope through a matrix zonotope,
following [AlthoffGK11; Proposition 4](@citet).

### Input

- `lm` -- a linear map of a zonotope through a matrix zonotope

### Output

A zonotope overapproximating the linear map.

"""
function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                    MAT<:MatrixZonotope{NM}}
    MZ = matrix(lm)
    Z = set(lm)
    T = promote_type(N, NM)

    n = dim(Z)
    w = ngens(MZ)
    h = ngens(Z)

    cZ = center(Z)
    GZ = genmat(Z)

    c = center(MZ) * cZ
    # generator
    G = Matrix{T}(undef, n, h * (w + 1) + w)
    G[:, 1:h] = center(MZ) * GZ
    @inbounds for (i, Ai) in enumerate(generators(MZ))
        G[:, (h * i + 1):(h * (i + 1))] = Ai * GZ
        G[:, h * (w + 1) + i] = Ai * cZ
    end

    return Zonotope(c, G)
end

"""
    overapproximate(lm::LinearMap{N,S,NM,MAT},
                         ::Type{<:Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                    MAT<:MatrixZonotopeProduct{NM}}

Overapproximate the linear map of a zonotope through a product of matrix zonotopes,
by recursively applying the overapproximation rule from the inside out.

### Input

- `lm` -- a linear map of a zonotope through a `MatrixZonotopeProduct`
- `Zonotope` -- target type

### Output

An overapproximation of the linear map as a zonotope.
"""
function overapproximate(lm::LinearMap{N,S,NM,MAT},
                         T::Type{<:Zonotope}) where {N,S<:AbstractZonotope{N},NM,
                                                     MAT<:MatrixZonotopeProduct{NM}}
    MZs = factors(matrix(lm))
    Z = set(lm)
    reduced = foldr((A, acc) -> overapproximate(A * acc, T), MZs; init=Z)
    return reduced
end
