using .Javis: Luxor, background, sethue, Video, Background, Object, render

export ground,
    luxify,
    nfolds

# default canvas background
function ground(args...)
    background("white") # canvas background
    return sethue("black") # pen color
end

# obtain a polygonal overapproximation with the given tolerance and then return
# a vector of `Luxor.Point`s
function luxify(X::LazySet; tol=1e-3)
    @assert dim(X) == 2 "can only work with two-dimensional sets"

    # ε-close outer polygonal approximation of X
    P = overapproximate(X, HPolygon, tol)

    # convert to vertex representation
    vlist = vertices_list(P)

    out_luxor = [Luxor.Point(Tuple(p)) for p in vlist]

    #return Luxor.poly(out_luxor, :stroke, close=true)
    return out_luxor
end

# useful function to split the array x into n equal parts (if possible)
# if the `rev` flag to return the partition in reverse order
# this function is useful to split in equal parts the animation time for an array of sets
function nfolds(x::AbstractVector, n::Int; rev=true)
    p = length(x)
    m = min(n, p)
    s = p / m
    out = [x[(round(Int64, (i - 1) * s) + 1):min(p, round(Int64, i * s))] for i in 1:m]
    return rev ? reverse(out) : out
end

# animate an array of sets using Javis.jl
# TODO pass colors, eg. with default cols = Javis.distinguishable_colors(length(X))
function Javis.animate(X::AbstractVector{ST}; kwargs...) where {N,ST<:LazySet{N}}

    # size of the animation
    sz = get(kwargs, :size, (500, 500))

    # number of sets to render
    nsets = length(X)

    # maximum time of the animation
    tmax = get(kwargs, :tmax, 10 * nsets)

    # number of time divisions (pass either an integer or a partition)
    ntdiv = get(kwargs, :ntdiv, nsets)
    tdiv = nfolds(1:tmax, ntdiv)

    # convert to Luxor.jl object
    Xlux = luxify.(X)

    # pick scaling factor
    α = get(kwargs, :scale, 50.0) # FIXME what is a good default? use radius of box_approximation(X[1]) ?

    # create empty video
    vid = Video(sz...)
    Background(1:tmax, ground)

    # add objects to the canvas
    for (dt, p) in zip(tdiv, Xlux)
        Object(dt, (args...) -> Javis.poly(α .* p, :stroke; close=true))
    end

    # render video
    return render(vid)
end
