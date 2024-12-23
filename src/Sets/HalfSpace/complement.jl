"""
# Extended help

    complement(H::HalfSpace)

### Algorithm

The result is the complementary half-space to `H`. If ``H: ⟨ a, x ⟩ ≤ b``,
then this method returns the half-space ``H′: ⟨ a, x ⟩ ≥ b``.
(Note that complementarity is understood in a relaxed sense, since the
intersection of ``H`` and ``H′`` is non-empty).
"""
function complement(H::HalfSpace)
    return HalfSpace(-H.a, -H.b)
end
