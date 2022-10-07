function writevtk(X::LazySet; file="output")
    points, connec = triangulate(X)
    ntriangles = length(connec)
    cell = VTKPolyhedron(1:ntriangles, connec...)

    vtk_grid(file, points, [cell]; compress=false) do vtk
        #
    end
end

function writevtk(X::AbstractVector{<:LazySet}; file="output")
    points_list = Vector{Matrix{Float32}}()
    cell_list = Vector{VTKPolyhedron}()
    count = 0
    for Xi in X
        points, connec = triangulate(Xi)
        ntriangles = length(connec)
        connec = [c .+ count for c in connec]
        cell = VTKPolyhedron(1:ntriangles, connec...)
        count += size(points, 2)
        push!(points_list, points)
        push!(cell_list, cell)
    end
    points_matrix = reduce(hcat, points_list)

    vtk_grid(file, points_matrix,  cell_list; compress=false) do vtk
        #
    end
end
