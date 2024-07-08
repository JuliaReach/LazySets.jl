"""
    project(S::Singleton, block::AbstractVector{Int}; kwargs...)

Concrete projection of a singleton.

### Input

- `S`     -- singleton
- `block` -- block structure, a vector with the dimensions of interest

### Output

A set representing the projection of the singleton `S` on the dimensions
specified by `block`.
"""
function project(S::Singleton, block::AbstractVector{Int}; kwargs...)
    return Singleton(element(S)[block])
end
