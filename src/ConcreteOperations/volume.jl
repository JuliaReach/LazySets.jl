export volume

"""
    volume(H::AbstractHyperrectangle)

Return the volume of a hyperrectangular set.

### Input

- `H` -- hyperrectangular set

### Output

The (Euclidean) volume of ``H``.

### Algorithm

The volume of the hyperrectangle ``H`` with vector radius ``r`` is
``2ⁿ ∏ᵢ rᵢ`` where ``rᵢ`` denotes the ``i``-th component of ``r``.
"""
function volume(H::AbstractHyperrectangle{N}) where {N<:AbstractFloat}
    r = radius_hyperrectangle(H)
    n = dim(H)
    α = exp(n*log(2))
    vol = α * prod(r)
    return vol
end

# fallback for rational
function volume(H::AbstractHyperrectangle{N}) where {N}
    r = radius_hyperrectangle(H)
    vol = prod(2 .* r)
    return vol
end

"""
    volume(X::LazySet{N}; partition=fill(2, dim(X)), tol=1e-2, itermax=5,
           inclusion_check_algorithm::String="constraints",
           disjointness_check_algorithm::String="exact") where {N}

Return an overapproximation of the volume of a set by partitioning into boxes.

### Input

- `X` -- set
- `partition` -- (optional, default: `fill(2, dim(X))`) scheme to partition the search
                 space at each iteration
- `tol`       -- (optional, default: `1e-2`) tolerance criterion for termination
- `itermax`   -- (optional, default: `5`) number of iterations criterion for termination
- `inclusion_check_algorithm` -- (optional, default: `"constraints"`) algorithm
                 for the inclusion check; valid options are `"constraints"` and
                 `"vertices"`
- `disjointness_check_algorithm` -- (optional, default: `"exact"`) algorithm for
                 the disjointness check; valid options are `"exact"` and `"sufficient"`

### Output

The volume of ``X``.

### Algorithm

TO-DO
"""
function volume(X::LazySet{N}; partition=fill(2, dim(X)), tol=1e-2, itermax=5,
                inclusion_check_algorithm::String="constraints",
                disjointness_check_algorithm::String="exact") where {N}

    # initial box overapproximation
    H = overapproximate(X, Hyperrectangle)

    # if X is box-shaped, its volume is that of H
    if issubset(H, X, algorithm=inclusion_check_algorithm)
        return volume(H)
    end

    # sequence of inner and outer volume computations
    vol_inner = Vector{N}()
    vol_outer = Vector{N}()

    # initial splitting
    Hsplit = split(H, partition)
    α = prod(partition)

    # iterative splitting
    k = 1
    while k <= itermax
        #println("iteration $k")
        #println("number of elements in the list $(length(Hsplit))")
        discard_and_evaluate!(Hsplit, vol_inner, vol_outer, X,
            inclusion_check_algorithm, disjointness_check_algorithm)
        bisect!(Hsplit, partition, α)
        if k >= 2 && (last(vol_outer) + vol_inner[end] - vol_inner[end-1] < tol)
            break
        end
        k += 1
    end

    return sum(vol_inner) + last(vol_outer)
end

# discard boxes which either are totally outside X, or totally inside X
# we also evaluate the volume for the inner boxes, and for the boxes whose intersection
# with X is non-empty but are not contained
function discard_and_evaluate!(Hsplit, vol_inner::Vector{N}, vol_outer::Vector{N}, X,
                               inclusion_check_algorithm,
                               disjointness_check_algorithm) where {N}
    vol_inner_current = zero(N)
    vol_outer_current = zero(N)
    delete_idx = Vector{Int}()
    for k in eachindex(Hsplit)
        Hk = Hsplit[k]
        if isdisjoint(Hk, X, algorithm=disjointness_check_algorithm)
            # discard this element
            push!(delete_idx, k)
        elseif issubset(H, X, algorithm=inclusion_check_algorithm)
            # accumulate volume and discard this element
            vol_inner_current += volume(Hk)
            push!(delete_idx, k)
        else
            # volume over-estimation by adding the volume of the boxes partially intersecting X
            vol_outer_current += volume(Hk)
        end
    end
    deleteat!(Hsplit, delete_idx)
    push!(vol_inner, vol_inner_current)
    push!(vol_outer, vol_outer_current)
end

# bisect all elements in the waiting list
function bisect!(Hsplit, partition, α)
    # length of Hsplit for the current iteration
    N = length(Hsplit)
    for i in 1:N
        μ = α * (i-1)
        Hi = split(Hsplit[1 + μ], partition)
        deleteat!(Hsplit, 1 + μ)
        for (k, x) in enumerate(Hi)
            insert!(Hsplit, k + μ, x)
        end
    end
end
