function Base.convert(::Type{<:VTKPolyhedron}, X::LazySet; shift=0)

    dim(X) == 3 || throw(ArgumentError("the dimension of the set should be three, got $(dim(X))"))

    poly = polyhedron(convert(HPolytope, X))
    mes = Mesh(poly)
    coords = Polyhedra.GeometryBasics.coordinates(mes)
    connec = Polyhedra.GeometryBasics.faces(mes)

    ntriangles = length(connec)
    npoints = 3*ntriangles
    points = Matrix{Float32}(undef, 3, npoints)

    for i in 1:npoints
        points[:, i] .= coords[i].data
    end

    connec_tup = getfield.(connec, :data)
    if !iszero(shift)
         connec_tup = [c .+ shift for c in connec_tup]
    end
    cell = VTKPolyhedron(1:ntriangles, connec_tup...)

    return points, cell
end

function writevtk(X::LazySet{N}; file="output") where {N}
    points, cell = convert(VTKPolyhedron, X)
    vtk_grid(file, points, [cell,]; compress = false) do vtk
        #
    end
end

function writevtk(X::AbstractVector{VT}; file="output") where {N, VT<:LazySet{N}}
    points_list = Vector{Matrix{Float32}}()
    cell_list = Vector{VTKPolyhedron}()
    count = 0
    for Xi in X
        points, cell = convert(VTKPolyhedron, Xi, shift=count)
        count += size(points, 2)
        push!(points_list, points)
        push!(cell_list, cell)
    end
    points_matrix = reduce(hcat, points_list)

    vtk_grid(file, points_matrix,  cell_list; compress = false) do vtk
        #
    end
end
