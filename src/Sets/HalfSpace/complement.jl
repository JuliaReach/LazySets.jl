"""
    complement(H::HalfSpace)

Return the complement of a half-space.

### Input

- `H` -- half-space

### Output

The half-space that is complementary to `H`. If ``H: ⟨ a, x ⟩ ≤ b``,
then this function returns the half-space ``H′: ⟨ a, x ⟩ ≥ b``.
(Note that complementarity is understood in a relaxed sense, since the
intersection of ``H`` and ``H′`` is non-empty).
"""
function complement(H::HalfSpace)
    return HalfSpace(-H.a, -H.b)
end
