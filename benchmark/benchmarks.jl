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

SUITE = BenchmarkGroup()  # parent BenchmarkGroup to contain our suite

include("support_vector.jl")
