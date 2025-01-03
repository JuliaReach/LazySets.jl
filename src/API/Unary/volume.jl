"""
    volume(X::LazySet)

Compute the volume, or Lebesgue measure, of a set.

### Input

- `X` -- set

### Output

A real number representing the Lebesgue measure of `X`.

### Notes

The [Lebesgue measure](https://en.wikipedia.org/wiki/Lebesgue_measure) has the
following common special cases:

- In 1D, it coincides with the *length*.
- In 2D, it coincides with the *area* (see also [`area`](@ref)).
- In 3D, it coincides with the *volume*.

In higher dimensions, it is also known as the *hypervolume* or simply *volume*.
"""
function volume(::LazySet) end
