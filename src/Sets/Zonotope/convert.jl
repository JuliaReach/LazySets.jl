"""
    convert(::Type{Zonotope}, Z::AbstractZonotope)

Convert a zonotopic set to a zonotope.

### Input

- `Zonotope` -- target type
- `H`        -- zonotopic set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, Z::AbstractZonotope)
    return Zonotope(center(Z), genmat(Z))
end
