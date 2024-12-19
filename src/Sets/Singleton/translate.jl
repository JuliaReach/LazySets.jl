"""
# Extended help

    translate!(S::Singleton, v::AbstractVector)

### Algorithm

We add the vector to the point in the singleton.
"""
function translate!(S::Singleton, v::AbstractVector)
    @assert length(v) == dim(S) "cannot translate a $(dim(S))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    S.element .+= v
    return S
end
