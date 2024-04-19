# convenience method for IntervalArithmetic.Interval
function vertices_list(x::IA.Interval{N}) where {N}
    a = IA.inf(x)
    b = IA.sup(x)
    ST = IA.SVector{1,N}
    return _isapprox(a, b) ? [ST(a)] : [ST(a), ST(b)]
end

# convenience method for IntervalArithmetic.IntervalBox
function vertices_list(H::IA.IntervalBox)
    return vertices_list(convert(Hyperrectangle, H))
end

const zero_itv(N) = IA.interval(zero(N))
const sym_itv(N) = IA.interval(-one(N), one(N))

zero_box(n::Int, N=Float64) = IA.IntervalBox(zero_itv(N), n)
sym_box(n::Int, N=Float64) = IA.IntervalBox(sym_itv(N), n)

"""
    fast_interval_pow(a::IntervalArithmetic.Interval, n::Int)

Compute the `n`th power of an interval without using correct rounding.

### Input

- `a` -- interval (from `IntervalArithmetic.jl`)
- `n` -- integer

### Output

A non-rigorous approximation of `a^n`.

### Notes

For a rigorous approximation with correct rounding, use `a^n` from
`IntervalArithmetic.jl`.
"""
function fast_interval_pow(a::IA.Interval, n::Int)
    # TODO review after IntervalArithmetic.jl#388
    if iszero(n)
        return one(a)
    elseif isodd(n)
        return IA.interval(a.lo^n, a.hi^n)
    else
        if 0 ∈ a
            return IA.interval(zero(a.lo), max(abs(a.lo), abs(a.hi))^n)
        else
            lon = a.lo^n
            hin = a.hi^n
            return IA.interval(min(lon, hin), max(lon, hin))
        end
    end
end
