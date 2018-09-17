using MathProgBase, GLPKMathProgInterface

export VPolytope,
       vertices_list

"""
    VPolytope{N<:Real} <: AbstractPolytope{N}

Type that represents a convex polytope in V-representation.

### Fields

- `vertices` -- the list of vertices
"""
struct VPolytope{N<:Real} <: AbstractPolytope{N}
    vertices::Vector{Vector{N}}
end

# constructor for a VPolytope with empty vertices list
VPolytope{N}() where {N<:Real} = VPolytope{N}(Vector{Vector{N}}())

# constructor for a VPolytope with no vertices of type Float64
VPolytope() = VPolytope{Float64}()

# constructor from a polygon in V-representation
function VPolytope(P::VPolygon, share::Bool=false)
    v = vertices_list(P)
    if !share
        v = copy(v)
    end
    return VPolytope(v)
end


# --- LazySet interface functions ---


"""
    dim(P::VPolytope)::Int

Return the dimension of a polytope in V-representation.

### Input

- `P`  -- polytope in V-representation

### Output

The ambient dimension of the polytope in V-representation.
If it is empty, the result is ``-1``.

### Examples

```jldoctest
julia> v = VPolytope();

julia> dim(v) > 0
false

julia> v = VPolytope([ones(3)])
VPolytope{Float64}(Array{Float64,1}[[1.0, 1.0, 1.0]])

julia> dim(v) == 3
true

```
"""
function dim(P::VPolytope)::Int
    return length(P.vertices) == 0 ? -1 : length(P.vertices[1])
end

"""
    σ(d::AbstractVector{N}, P::VPolytope{N}; algorithm="hrep") where {N<:Real}

Return the support vector of a polyhedron (in V-representation) in a given
direction.

### Input

- `d`         -- direction
- `P`         -- polyhedron in V-representation
- `algorithm` -- (optional, default: `'hrep'`) method to compute the support vector

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N},
           P::VPolytope{N};
           algorithm="hrep") where {N<:Real}
    if algorithm == "hrep"
        @assert isdefined(@__MODULE__, :Polyhedra) "this algorithm needs the " *
            "package 'Polyhedra' loaded"
        return σ(d, tohrep(P))
    else
        error("the algorithm $algorithm is not known")
    end
end

# --- VPolytope functions ---

"""
    vertices_list(P::VPolytope{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polytope in V-representation.

### Input

- `P` -- polytope in vertex representation

### Output

List of vertices.
"""
function vertices_list(P::VPolytope{N})::Vector{Vector{N}} where {N<:Real}
    return P.vertices
end


function load_polyhedra_vpolytope() # function to be loaded by Requires
return quote
# see the interface file AbstractPolytope.jl for the imports

export intersection, convex_hull, cartesian_product, vertices_list, tohrep, tovrep

# VPolytope from a VRep
function VPolytope(P::VRep{N, T}, backend=default_polyhedra_backend(N)) where {N, T}
    vertices = Vector{Vector{T}}()
    for vi in Polyhedra.points(P)
        push!(vertices, vi)
    end
    return VPolytope(vertices)
end

"""
    polyhedron(P::VPolytope{N}, [backend]=default_polyhedra_backend(N)) where {N}

Return an `VRep` polyhedron from `Polyhedra.jl` given a polytope in V-representation.

### Input

- `P`       -- polytope
- `backend` -- (optional, default: `default_polyhedra_backend(N)`) the polyhedral
               computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
               for further information

### Output

A `VRep` polyhedron.
"""
function polyhedron(P::VPolytope{N}, backend=default_polyhedra_backend(N)) where {N}
    V = hcat(vertices_list(P)...)'
    return polyhedron(Polyhedra.vrep(V), backend)
end

"""
    intersection(P1::VPolytope{N}, P2::VPolytope{N};
                [backend]=default_polyhedra_backend(N),
                [prunefunc]=removehredundancy!)::VPolytope{N} where {N<:Real}

Compute the intersection of two polytopes in V-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`) the polyhedral
                  computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information
- `prunefunc` -- (optional, default: `removehredundancy!`) function to post-process
                  the output of `intersect`

### Output

The `VPolytope` obtained by the intersection of `P1` and `P2`.
"""
function intersection(P1::VPolytope{N}, P2::VPolytope{N};
                      backend=default_polyhedra_backend(N),
                      prunefunc=removehredundancy!)::VPolytope{N} where {N<:Real}

    P1 = polyhedron(P1, backend)
    P2 = polyhedron(P2, backend)
    Pint = Polyhedra.intersect(P1, P2)
    prunefunc(Pint)
    return VPolytope(Pint)
end

"""
    convex_hull(P1::VPolytope{N}, P2::VPolytope{N}; backend=default_polyhedra_backend(N)) where {N}

Compute the convex hull of the set union of two polytopes in V-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`) the polyhedral
                  computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `VPolytope` obtained by the concrete convex hull of `P1` and `P2`.
"""
function convex_hull(P1::VPolytope{N}, P2::VPolytope{N}; backend=default_polyhedra_backend(N)) where {N}
    Pch = convexhull(polyhedron(P1, backend), polyhedron(P2, backend))
    return VPolytope(Pch)
end

"""
    cartesian_product(P1::VPolytope{N}, P2::VPolytope{N}; backend=default_polyhedra_backend(N)) where {N}

Compute the Cartesian product of two polytopes in V-representation.

### Input

- `P1`         -- polytope
- `P2`         -- another polytope
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`) the polyhedral
                  computations backend, see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `VPolytope` obtained by the concrete Cartesian product of `P1` and `P2`.
"""
function cartesian_product(P1::VPolytope{N}, P2::VPolytope{N}; backend=default_polyhedra_backend(N)) where {N}
    Pcp = hcartesianproduct(polyhedron(P1, backend), polyhedron(P2, backend))
    return VPolytope(Pcp)
end

"""
    tohrep(P::VPolytope{N}; backend=default_polyhedra_backend(N)) where {N}

Transform a polytope in V-representation to a polytope in H-representation.

### Input

- `P`          -- polytope in vertex representation
- `backend`    -- (optional, default: `default_polyhedra_backend(N)`) the polyhedral
                  computations backend,
                  see [Polyhedra's documentation](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)
                  for further information

### Output

The `HPolytope` which is the constraint representation of the given polytope
in vertex representation.
"""
function tohrep(P::VPolytope{N}; backend=default_polyhedra_backend(N)) where {N}
    P = polyhedron(P, backend)
    return HPolytope(P)
end

"""
    tovrep(P::VPolytope)

Return a vertex representation of the given polytope in vertex
representation (no-op).

### Input

- `P` -- polytope in vertex representation

### Output

The same polytope instance.
"""
function tovrep(P::VPolytope)
    return P
end

end # quote
end # function load_polyhedra_vpolytope()
