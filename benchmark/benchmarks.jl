#=
This file defines a benchmark suite with the tools provided by
`PkgBenchmark` and `BenchmarkTools`.

To run the benchmarks, execute:

```julia
julia> using PkgBenchmark
julia> results = benchmarkpkg("LazySets")
```

To compare current version to another tagged version, commit or branch:

```julia
julia> results = judge("LazySets", <tagged-version-or-branch>)
```

To export the benchmark results to a Markwodn file:

```julia
julia> export_markdown("results.md", results)
```

To export the benchmark results to a JSON file:

```julia
julia> writeresults("results.json", results)
```
=#
using BenchmarkTools, LazySets

function tridiagm(a, b, c, n)  # to generate tridiagonal matrices
       dd, du = ones(n), ones(n - 1)
       b*diagm(dd) + a*diagm(du, -1) + c*diagm(du, 1)
end

SUITE = BenchmarkGroup()  # parent BenchmarkGroup to contain our suite

SUITE["Balls"] = BenchmarkGroup()
SUITE["Linear map"] = BenchmarkGroup()
begin
    for set_type in (Ball1, Ball2, BallInf)
        for n in (1, 10, 100, 1000)
            X = set_type(ones(n), 0.1)
            d = [i/n for i in 1:n]
            SUITE["Balls"][string(set_type), n] = @benchmarkable σ($d, $X)

            A = tridiagm(1, -2, 0.5, n)
            Y = A * X
            SUITE["Linear map"][string(set_type), n] = @benchmarkable σ($d, $Y)
        end
    end
end
