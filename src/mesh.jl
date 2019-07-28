function load_mesh()
return quote

using .Polyhedra: Mesh
using .Makie: mesh, mesh!
import .Makie.AbstractPlotting: Automatic

export plot3d, plot3d!

function plot3d_helper(S::LazySet{N}, backend) where {N}
    @assert dim(S) <= 3 "plot3d can only be used to plot sets of dimension three (or lower); " *
        "but the given set is $(dim(S))-dimensional"

    @assert applicable(constraints_list, S) "plot3d requires that the list of constraints of `S`, " *
        "`constraints_list(S)` is applicable; try overapproximating with an `HPolytope` first"

    P = HPolytope(constraints_list(S))
    remove_redundant_constraints!(P)
    P_poly = polyhedron(P, backend=backend)
    P_poly_mesh = Mesh(P_poly)
    return P_poly_mesh
end

function plot3d(S::LazySet{N}; backend=default_polyhedra_backend(S, N),
                alpha=1.0, color=:blue, colormap=:viridis, colorrange=Automatic(), interpolate=false,
                linewidth=1, overdraw=false, shading=true, transparency=true, visible=true) where {N}

    P_poly_mesh = plot3d_helper(S, backend)
    return mesh(P_poly_mesh, alpha=alpha, color=color, colormap=colormap, colorrange=colorrange,
                interpolate=interpolate, linewidth=linewidth, transparency=transparency, visible=visible)
end

function plot3d!(S::LazySet{N}; backend=default_polyhedra_backend(S, N),
                alpha=1.0, color=:blue, colormap=:viridis, colorrange=Automatic(), interpolate=false,
                linewidth=1, overdraw=false, shading=true, transparency=true, visible=true) where {N}

    P_poly_mesh = plot3d_helper(S, backend)
    return mesh!(P_poly_mesh, alpha=alpha, color=color, colormap=colormap, colorrange=colorrange,
                 interpolate=interpolate, linewidth=linewidth, transparency=transparency, visible=visible)
end

end # quote
end # function load_mesh()
