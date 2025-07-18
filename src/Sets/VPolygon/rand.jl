# range of the default number of vertices in `rand`
const DEFAULT_RAND_VERTEX_RANGE = 3:10

"""
# Extended help

    rand(::Type{VPolygon}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Input

- `num_vertices` -- (optional, default: `-1`) number of vertices of the
                    polygon (see comment below)

### Notes

The number of vertices can be controlled with the argument `num_vertices`.
For a negative value we choose a random number in the range
`$DEFAULT_RAND_VERTEX_RANGE`.

### Algorithm

We follow the idea described [here](https://stackoverflow.com/a/47358689) based
on [Valtr95](@citet). There is also a nice video available
[here](http://cglab.ca/~sander/misc/ConvexGeneration/convex.html).
"""
function rand(::Type{VPolygon};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              num_vertices::Int=-1)
    @assert dim == 2 "cannot create a random VPolygon of dimension $dim"
    rng = reseed!(rng, seed)
    if num_vertices < 0
        num_vertices = rand(rng, DEFAULT_RAND_VERTEX_RANGE)
    end

    # special cases, 0 or 1 vertex
    if num_vertices == 0
        return VPolygon{N}()
    elseif num_vertices == 1
        return VPolygon([randn(rng, N, 2)])
    end

    # general case, >= 2 vertices

    # get random horizontal and vertical vectors
    horiz = rand_pos_neg_zerosum_vector(num_vertices; N=N, rng=rng)
    vert = rand_pos_neg_zerosum_vector(num_vertices; N=N, rng=rng)

    # randomly combine horizontal and vertical vectors
    m = num_vertices
    directions = Vector{Vector{N}}(undef, num_vertices)
    shuffle(rng, vert)
    for (i, x) in enumerate(horiz)
        y = splice!(vert, rand(rng, 1:m))
        directions[i] = [x, y]
        m -= 1
    end
    sort!(directions; lt=<=) # sort by angle

    # connect directions
    vertices = Vector{Vector{N}}(undef, num_vertices)
    # random starting point
    @inbounds begin
        vertices[1] = randn(rng, N, 2)
        for i in 1:(length(directions) - 1)
            vertices[i + 1] = vertices[i] + directions[i]
        end
        @assert _isapprox(vertices[end] + directions[end], vertices[1])
    end
    return VPolygon(vertices; apply_convex_hull=true)
end
