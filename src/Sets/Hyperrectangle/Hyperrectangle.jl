"""
    Hyperrectangle{N, VNC<:AbstractVector{N}, VNR<:AbstractVector{N}
                  } <: AbstractHyperrectangle{N}

Type that represents a hyperrectangle.

A [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) is the
Cartesian product of one-dimensional intervals.

### Fields

- `center` -- center of the hyperrectangle as a real vector
- `radius` -- radius of the hyperrectangle as a real vector, i.e., half of its
              width along each coordinate direction

### Examples

The `Hyperrectangle` type stores a vector representing the center and another
vector representing the radius. The default constructor `Hyperrectangle(c, r)`
receives the center and radius, in that order. For instance,

```jldoctest hyperrectangle_constructor
julia> c = [-1.0, 1.0];

julia> r = [2.0, 1.0];

julia> H = Hyperrectangle(c, r)
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([-1.0, 1.0], [2.0, 1.0])
```

The above instance represents the hyperrectangle with the following vertices:

```jldoctest hyperrectangle_constructor
julia> vertices_list(H)
4-element Vector{Vector{Float64}}:
 [1.0, 2.0]
 [-3.0, 2.0]
 [1.0, 0.0]
 [-3.0, 0.0]
```

The getter functions for the center and the radius are `center` and
`radius_hyperrectangle` (since `radius` corresponds to the radius of the
enclosing ball of minimal volume):

```jldoctest hyperrectangle_constructor
julia> center(H)
2-element Vector{Float64}:
 -1.0
  1.0

julia> radius_hyperrectangle(H)
2-element Vector{Float64}:
 2.0
 1.0
```

There is also a constructor from lower and upper bounds with keyword arguments
`high` and `low`. The following construction results in the same hyperrectangle
as in the previous paragraph:

```jldoctest hyperrectangle_constructor
julia> l = [-3.0, 0.0];

julia> h = [1.0, 2.0];

julia> Hyperrectangle(low=l, high=h)
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([-1.0, 1.0], [2.0, 1.0])
```

By default, the constructor checks that that radius of the hyperrectangle
is nonnegative. To suppress this check, use the `check_bounds` optional flag
in the constructor. Note that if `check_bounds` is set to `false`, the behavior
of a set with contradictory bounds is undefined.
"""
struct Hyperrectangle{N,VNC<:AbstractVector{N},VNR<:AbstractVector{N}} <: AbstractHyperrectangle{N}
    center::VNC
    radius::VNR

    # default constructor with length comparison & domain constraint for radius
    function Hyperrectangle(center::VNC, radius::VNR;
                            check_bounds::Bool=true) where
             {N,VNC<:AbstractVector{N},VNR<:AbstractVector{N}}
        @assert length(center) == length(radius) "length of center and radius must be equal"
        if check_bounds
            @assert all(v -> v >= zero(N), radius) "radius must be nonnegative but is $radius"
        end
        return new{N,VNC,VNR}(center, radius)
    end
end

# constructor from keyword arguments (lower and upper bounds)
function Hyperrectangle(; high::AbstractVector,
                        low::AbstractVector,
                        check_bounds::Bool=true)
    # compute center and radius from high and low vectors
    center = (high .+ low) ./ 2
    radius = high .- center
    return Hyperrectangle(center, radius; check_bounds=check_bounds)
end
