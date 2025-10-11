function writevtk(X::LazySet; file="output")
    points, connection = triangulate_faces(X)
    ntriangles = length(connection)
    cell = VTKPolyhedron(1:ntriangles, connection...)

    vtk_grid(file, points, [cell]; compress=false) do vtk
        #
    end
end

function writevtk(X::AbstractVector{<:LazySet}; file="output")
    points_list = Vector{Matrix{Float32}}()
    cell_list = Vector{VTKPolyhedron}()
    count = 0
    for Xi in X
        points, connection = triangulate_faces(Xi)
        ntriangles = length(connection)
        connection = [c .+ count for c in connection]
        cell = VTKPolyhedron(1:ntriangles, connection...)
        count += size(points, 2)
        push!(points_list, points)
        push!(cell_list, cell)
    end
    points_matrix = reduce(hcat, points_list)

    vtk_grid(file, points_matrix, cell_list; compress=false) do vtk
        #
    end
end
