"""
# Extended help

    translate!(S::Singleton, v::AbstractVector)

### Algorithm

We add the vector to the point in the singleton.
"""
@validate function translate!(S::Singleton, v::AbstractVector)
    S.element .+= v
    return S
end
