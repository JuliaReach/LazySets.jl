eval(quote
using .Javis: Luxor, background, sethue

export ground,
       luxify,
       nfolds,
       animate

# default canvas background
function ground(args...)
    background("white") # canvas background
    sethue("black") # pen color
end

# convert a lazysets to a polygonal overapproximation with the given tolerance,
# and then to the function `poly` from Luxor.jl that receives array of points
function luxify(X::LazySet; tol=1e-3)

    # ε-close outer polygonal approximation of X
    P = overapproximate(X, HPolygon, tol)

    # convert to vertex representation
    V = convert(VPolygon, P)
    out = vertices_list(P)

    out_luxor = [Luxor.Point(Tuple(p)) for p in out]

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
    out = [x[round(Int64, (i-1)*s)+1:min(p,round(Int64, i*s))] for i in 1:m]
    return rev ? reverse(out) : out
end

# animate an array of sets using Javis.jl
# TODO pass colors, eg. with default cols = Javis.distinguishable_colors(length(X))
function animate(X::AbstractVector{ST}; kwargs...) where {N, ST<:LazySet{N}}

    # size of the animation
    sz = get(kwargs, :size, (500, 500))

    # number of sets to render
    nsets = length(X)

    # maximum time of the animation
    tmax = get(kwargs, :tmax, 10*nsets)

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
        Object(dt, (args...) -> poly(α .* p, :stroke, close=true))
    end

    # render video
    render(vid)
end

end)
