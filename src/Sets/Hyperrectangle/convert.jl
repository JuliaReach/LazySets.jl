"""
    convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)

Convert a hyperrectangular set to a hyperrectangle.

### Input

- `Hyperrectangle` -- hyperrectangle target type
- `H`              -- hyperrectangular set

### Output

A hyperrectangle.

### Examples

```jldoctest
julia> convert(Hyperrectangle, Interval(0.0, 1.0))
Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}([0.5], [0.5])
```
"""
function convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)
    return Hyperrectangle(center(H), radius_hyperrectangle(H))
end

function load_IntervalArithmetic_convert()
    return quote
        import IntervalArithmetic as IA

        """
            convert(::Type{IntervalArithmetic.IntervalBox}, H::AbstractHyperrectangle)

        Convert a hyperrectangular set to an `IntervalBox` from `IntervalArithmetic`.

        ### Input

        - `IntervalBox` -- target type
        - `H`           -- hyperrectangular set

        ### Output

        An `IntervalBox`.
        """
        function convert(::Type{IA.IntervalBox}, H::AbstractHyperrectangle)
            return IA.IntervalBox(IA.interval.(low(H), high(H)))
        end

        """
            convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)

        Convert an `IntervalBox` from `IntervalArithmetic` to a hyperrectangular set.

        ### Input

        - `Hyperrectangle` -- target type
        - `IB`             -- interval box

        ### Output

        A `Hyperrectangle`.

        ### Notes

        `IntervalArithmetic.IntervalBox` uses *static* vectors to store each component
        interval; hence the resulting `Hyperrectangle` has its center and radius
        represented as a static vector (`SArray`).
        """
        function convert(::Type{Hyperrectangle}, IB::IA.IntervalBox)
            low_IB = IA.inf.(IB)
            high_IB = IA.sup.(IB)
            return Hyperrectangle(; low=low_IB, high=high_IB)
        end

        # method for Interval
        function convert(::Type{Hyperrectangle}, I::IA.Interval)
            low_I = [IA.inf(I)]
            high_I = [IA.sup(I)]
            return Hyperrectangle(; low=low_I, high=high_I)
        end
    end
end  # quote / load_IntervalArithmetic_convert
