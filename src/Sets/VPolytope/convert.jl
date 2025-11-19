"""
    convert(::Type{VPolytope}, X::LazySet; [prune]::Bool=true)

Convert a polytopic set to a polytope in vertex representation.

### Input

- `VPolytope` -- target type
- `X`         -- polytopic set
- `prune`     -- (optional, default: `true`) option to remove redundant vertices

### Output

The given set represented as a polytope in vertex representation.

### Algorithm

This method uses `vertices_list`. Use the option `prune` to select whether to
remove redundant vertices before constructing the polytope.
"""
function convert(::Type{VPolytope}, X::LazySet; prune::Bool=true)
    # TODO there should be a better mechanism to choose the right method
    if hasmethod(vertices_list, tuple(typeof(X)), (:prune,))
        vlist = vertices_list(X; prune=prune)
    else
        vlist = vertices_list(X)
    end
    return VPolytope(vlist)
end

function load_Polyhedra_convert_VPolytope()
    return quote
        using .Polyhedra: VRep

        # VPolytope from a VRep
        function convert(::Type{VPolytope}, P::VRep)
            vertices = collect(Polyhedra.points(P))
            return VPolytope(vertices)
        end
    end
end  # load_Polyhedra_convert_VPolytope
