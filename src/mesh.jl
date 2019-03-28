#using .LazySets
#using .LazySets.Approximations
#using LazySets: default_polyhedra_backend, dim # needed?

using Makie: mesh, mesh!
using AbstractPlotting: Automatic
using Polyhedra: Mesh

export plot3d
#import Plots: plot3d

function plot3d(S::LazySet{N}, num_dirs=10; backend=default_polyhedra_backend(S, N),
                alpha=1.0, color=:blue, colormap=:viridis, colorrange=Automatic(), interpolate=false,
                linewidth=1, overdraw=false, shading=true, transparency=true, visible=true) where {N}

    if dim(S) != 3
        throw(ArgumentError("plot3d can only be used to plot 3D sets, but the given set is $(dim(S))-dimensional"))
    end
    
    if !applicable(constraints_list, S)
        throw(ArgumentError("the `constraints_list` function is not applicable to a set of type $S, hence plot3d cannot be used"))
    end

    P = overapproximate(S, SphericalDirections(num_dirs))
    remove_redundant_constraints!(P)
    return plot3d(P, backend=backend, alpha=alpha, color=color, colormap=colormap, colorrange=colorrange,
                     interpolate=interpolate, linewidth=linewidth, transparency=transparency, visible=visible)
end

function plot3d(S::AbstractPolytope{N}; backend=default_polyhedra_backend(S, N),
                alpha=1.0, color=:blue, colormap=:viridis, colorrange=Automatic(), interpolate=false,
                linewidth=1, overdraw=false, shading=true, transparency=true, visible=true) where {N}

    if dim(S) != 3
        throw(ArgumentError("plot3d can only be used to plot 3D sets, but the given set is $(dim(S))-dimensional"))
    end

    P = convert(HPolytope, S)
    P_poly = polyhedron(P, backend=backend)
    P_poly_mesh = Mesh(P_poly)

    mesh(P_poly_mesh, alpha=alpha, color=color, colormap=colormap, colorrange=colorrange,
                      interpolate=interpolate, linewidth=linewidth, transparency=transparency, visible=visible)
end

function plot3d!(S::LazySet{N}, num_dirs=10; backend=default_polyhedra_backend(S, N),
                alpha=1.0, color=:blue, colormap=:viridis, colorrange=Automatic(), interpolate=false,
                linewidth=1, overdraw=false, shading=true, transparency=true, visible=true) where {N}

    if dim(S) != 3
        throw(ArgumentError("plot3d can only be used to plot 3D sets, but the given set is $(dim(S))-dimensional"))
    end
    
    if !applicable(constraints_list, S)
        throw(ArgumentError("the `constraints_list` function is not applicable to a set of type $S, hence plot3d cannot be used"))
    end

    P = overapproximate(S, SphericalDirections(num_dirs))
    remove_redundant_constraints!(P)
    return plot3d!(P, backend=backend, alpha=alpha, color=color, colormap=colormap, colorrange=colorrange,
                     interpolate=interpolate, linewidth=linewidth, transparency=transparency, visible=visible)
end

function plot3d!(S::AbstractPolytope{N}; backend=default_polyhedra_backend(S, N),
                alpha=1.0, color=:blue, colormap=:viridis, colorrange=Automatic(), interpolate=false,
                linewidth=1, overdraw=false, shading=true, transparency=true, visible=true) where {N}

    if dim(S) != 3
        throw(ArgumentError("plot3d can only be used to plot 3D sets, but the given set is $(dim(S))-dimensional"))
    end

    P = convert(HPolytope, S)
    P_poly = polyhedron(P, backend=backend)
    P_poly_mesh = Mesh(P_poly)

    mesh!(P_poly_mesh, alpha=alpha, color=color, colormap=colormap, colorrange=colorrange,
                      interpolate=interpolate, linewidth=linewidth, transparency=transparency, visible=visible)
end
