export plot3d, plot3d!

function load_polyhedra_mesh()
return quote

using .Polyhedra: Mesh

end end  # quote / function load_polyhedra_mesh()

function load_makie()
return quote

using .Makie: mesh, mesh!
using .Makie: Automatic

end end  # quote / function load_makie()

# helper function for 3D plotting; converts S to a polytope in H-representation
function _plot3d_helper(S::LazySet, backend)
    @assert dim(S) <= 3 "plot3d can only be used to plot sets of dimension 3 " *
        "(or lower), but the given set is $(dim(S))-dimensional"

    @assert applicable(constraints_list, S) "plot3d requires that the list " *
        "of constraints of `S`, `constraints_list(S)` is applicable; try " *
        "overapproximating with an `HPolytope` first"

    @assert isbounded(S) "plot3d requires a bounded set"

    P = HPolytope(constraints_list(S))
    remove_redundant_constraints!(P)
    P_poly = polyhedron(P, backend=backend)
    P_poly_mesh = Mesh(P_poly)
    return P_poly_mesh
end

"""
    plot3d(S::LazySet; [backend]=default_polyhedra_backend(S), [alpha]=1.0,
           [color]=:blue, [colormap]=:viridis, [colorrange]=nothing,
           [interpolate]=false, [linewidth]=1, [overdraw]=false, [shading]=true,
           [transparency]=true, [visible]=true)

Plot a three-dimensional set using `Makie`.

### Input

- `S`            -- set
- `backend`      -- (optional, default: `default_polyhedra_backend(S)`) backend
                    for polyhedral computations
- `alpha`        -- (optional, default: `1.0`) float in `[0,1]`; the alpha or
                    transparency value
- `color`        -- (optional, default: `:blue`) `Symbol` or `Colorant`; the
                    color of the main plot element (markers, lines, etc.), which
                    can be a color symbol/string like `:red`
- `colormap`     -- (optional, default: `:viridis`) the color map of the main
                    plot; use `available_gradients()` to see which gradients are
                    available, which can also be used as `[:red, :black]`
- `colorrange`   -- (optional, default: `nothing`, which falls back to
                    `Makie.Automatic()`) a tuple `(min, max)` where `min` and
                    `max` specify the data range to be used for indexing the
                    `colormap`
- `interpolate`  -- (optional, default: `false`) a boolean for heatmap and
                    images; toggles color interpolation between nearby pixels
- `linewidth`    -- (optional, default: `1`) a number that specifies the width
                    of the line in `line` and `linesegments` plots
- `overdraw`     -- (optional, default: `false`)
- `shading`      -- (optional, default: `true`) a boolean that toggles shading
                    (for meshes)
- `transparency` -- (optional, default: `true`) if `true`, the set is
                    transparent, otherwise it is displayed as a solid object
- `visible`      -- (optional, default: `true`) a boolean that toggles
                    visibility of the plot

For a complete list of attributes and usage see
[Makie's documentation](http://makie.juliaplots.org/stable/plot-attributes).

### Notes

This plot recipe works by computing the list of constraints of `S` and converting
to a polytope in H-representation. Then, this polytope is transformed with
`Polyhedra.Mesh` and plotted using the `mesh` function.

If the function `constraints_list` is not applicable to your set `S`, try
overapproximation first; e.g. via

```julia
julia> Sapprox = overapproximate(S, SphericalDirections(10))

julia> using Polyhedra, GLMakie

julia> plot3d(Sapprox)
```
The number `10` above corresponds to the number of directions considered; for
better resolution use higher values (but it will take longer).

For efficiency consider using the `CDDLib` backend, as in

```julia
julia> using CDDLib

julia> plot3d(Sapprox, backend=CDDLib.Library())
```

### Examples

The functionality requires *both* `Polyhedra` and a `Makie` backend. After
loading `LazySets`, do `using Polyhedra, GLMakie` (or another Makie backend).

```julia
julia> using LazySets, Polyhedra, GLMakie

julia> plot3d(10 * rand(Hyperrectangle, dim=3))

julia> plot3d!(10 * rand(Hyperrectangle, dim=3), color=:red)
```
"""
function plot3d(S::LazySet; backend=default_polyhedra_backend(S), alpha=1.0,
                color=:blue, colormap=:viridis, colorrange=nothing,
                interpolate=false, linewidth=1, overdraw=false, shading=true,
                transparency=true, visible=true)
    require(@__MODULE__, [:Makie, :Polyhedra]; fun_name="plot3d")

    if colorrange == nothing
        colorrange = Automatic()
    end
    P_poly_mesh = _plot3d_helper(S, backend)
    return mesh(P_poly_mesh, alpha=alpha, color=color, colormap=colormap,
                colorrange=colorrange, interpolate=interpolate,
                linewidth=linewidth, transparency=transparency, visible=visible)
end

"""
    plot3d!(S::LazySet; backend=default_polyhedra_backend(S), [alpha]=1.0,
           [color]=:blue, [colormap]=:viridis, [colorrange]=nothing,
           [interpolate]=false, [linewidth]=1, [overdraw]=false, [shading]=true,
           [transparency]=true, [visible]=true)

Plot a three-dimensional set using Makie.

### Input

See `plot3d` for the description of the inputs. For a complete list of
attributes and usage see [Makie's
documentation](http://makie.juliaplots.org/stable/plot-attributes).

### Notes

See the documentation of `plot3d` for examples.
"""
function plot3d!(S::LazySet; backend=default_polyhedra_backend(S), alpha=1.0,
                 color=:blue, colormap=:viridis, colorrange=nothing,
                 interpolate=false, linewidth=1, overdraw=false, shading=true,
                 transparency=true, visible=true)
    require(@__MODULE__, [:Makie, :Polyhedra]; fun_name="plot3d!")

    if colorrange == nothing
        colorrange = Automatic()
    end
    P_poly_mesh = _plot3d_helper(S, backend)
    return mesh!(P_poly_mesh, alpha=alpha, color=color, colormap=colormap,
                 colorrange=colorrange, interpolate=interpolate,
                 linewidth=linewidth, transparency=transparency,
                 visible=visible)
end
