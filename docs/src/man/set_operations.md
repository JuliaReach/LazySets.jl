# Operations on sets

In this section we show which typical set operations this library supports.

```@contents
Pages = ["set_operations.md"]
Depth = 3
```

## Unary operations

The following table lists all operations (apart from `dim` and ``, which must be
supported by all types) that take one convex set as argument in the columns.
In the rows we list all set types, both the interfaces (where we abbreviate the
`Abstract` prefix), the basic set types, and the lazy set operations, each
sorted alphabetically.
The table entries have the following meaning.
- An "x" indicates that the operation is implemented for the respective set.
- "(⋅)" indicates that the operation is inherited from a supertype.


| type ↓ \ operation →         | radius | diameter | norm | an_element | ∈ |
|------------------------------|--------|----------|------|------------|---|
| **Interfaces**               |        |          |      |            |   |
| `LazySet`                    |        |          |      |            |   |
| `AHPolygon`                  |        |          |      |     x      | x |
| `AHyperrectangle`            |   x    |    x     |  x   |    (x)     | x |
| `APointSymmetric`            |        |          |      |     x      |   |
| `APointSymmetricPolytope`    |        |          |      |     x      |   |
| `APolygon`                   |        |          |      |            |   |
| `APolytope`                  |        |          |      |            |   |
| `ASingleton`                 |  (x)   |   (x)    | (x)  |     x      | x |
|                              |        |          |      |            |   |
| **Basic set types**          |        |          |      |            |   |
| `Ball1`                      |        |          |      |    (x)     | x |
| `Ball2`                      |        |          |      |    (x)     | x |
| `BallInf`                    |   x    |   (x)    | (x)  |    (x)     |(x)|
| `Ballp`                      |        |          |      |    (x)     | x |
| `Ellipsoid`                  |        |          |      |    (x)     | x |
| `EmptySet`                   |        |          |      |     x      | x |
| `HPolygon`/`HPolygonOpt`     |        |          |      |    (x)     |(x)|
| `HPolytope`                  |        |          |      |            |   |
| `Hyperrectangle`             |  (x)   |   (x)    | (x)  |    (x)     |(x)|
| `Singleton`                  |  (x)   |   (x)    | (x)  |    (x)     |(x)|
| `VPolygon`                   |        |          |      |     x      | x |
| `ZeroSet`                    |  (x)   |   (x)    | (x)  |    (x)     | x |
| `Zonotope`                   |        |          |      |    (x)     | x |
|                              |        |          |      |            |   |
| **Lazy set operation types** |        |          |      |            |   |
| `CartesianProduct`           |        |          |      |            | x |
| `CartesianProductArray`      |        |          |      |            | x |
| `ConvexHull`                 |        |          |      |            |   |
| `ConvexHullArray`            |        |          |      |            |   |
| `ExponentialMap`             |        |          |      |            | x |
| `ExponentialProjectionMap`   |        |          |      |            |   |
| `HalfSpace`                  |        |          |      |     x      | x |
| `Hyperplane`                 |        |          |      |     x      | x |
| `Intersection`               |        |          |      |            |   |
| `LinearMap`                  |        |          |      |            | x |
| `MinkowskiSum`               |        |          |      |            |   |
| `MinkowskiSumArray`          |        |          |      |            |   |
| `SymmetricIntervalHull`      |  (x)   |   (x)    | (x)  |    (x)     |(x)|


### `radius`

This function returns the radius of a set.
It is defined as the enclosing ball (of the given norm) of minimal volume.

### `diameter`

This function returns the diameter of a set.
It is defined as the enclosing ball (of the given norm) of minimal volume.
The implementation is inherited for all set types if the norm is the infinity
norm, in which case the result is defined as twice the radius.

### `norm`

This function returns the norm of a set.
It is defined as the norm of the enclosing ball (of the given norm) of minimal
volume.

### `an_element`

This function returns some element in the set.
Consecutive calls to this function typically return the same element.

### `∈`

This function checks containment of a given vector in the set.
The operator can be used in infix notation (`v ∈ S`) and in inverse operand
order (`S ∋ v`).
Alias: `in`


## Binary operations

The following table lists all operations that take two convex set as argument in
the entries.
In the rows we list all set types, both the interfaces (where we abbreviate the
`Abstract` prefix), the basic set types, and the lazy set operations, each
sorted alphabetically.
In the columns we also list the operations, but abbreviated.
The table entries consist of subsets of the following list of operations.
- "⊆" stands for the subset check operation.
- "∩=∅" stands for the intersection emptiness check operation.
- "(⋅)" indicates that the operation is inherited from a supertype.
- "[⋅]" indicates that a possible type ambiguity can occur but is resolved.
- "{⋅}" indicates that a more efficient implementation is used instead of the inherited one.
- "**!⋅!**" (bold face) indicates that a type ambiguity may occur.


| type ↓ \ type →              |LazyS|AHPon|AHrec|APSym|APSPol|APgon|APtop|ASingle| B1  | B2  |BInf | Bp  |Ellip|EmpSt|HPgon|HPtop|Hrect |Single |VPgon|ZeroSet|Zonot| CP  | CPA | CH  | CHA |EMap | EPM |HalfS|HypPl|Inter|LMap |MSum | MSA | SIH |
|------------------------------|-----|-----|-----|-----|------|-----|-----|-------|-----|-----|-----|-----|-----|-----|-----|-----|------|-------|-----|-------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| **Interfaces**               |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `LazySet`                    |     |     |  ⊆  |     |      |     |     |(⊆),∩=∅|     |     | (⊆) |     |     |     |     |     | (⊆)  | (⊆)   |     | (⊆)   |     |     |     |     |     |     |     |     |     |     |     |     |     | (⊆) |
| `AHPolygon`                  | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `AHyperrectangle`            | (⊆) | (⊆) |{⊆},∩=∅|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |(⊆,∩=∅)|(⊆,∩=∅)|(⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `APointSymmetric`            |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `APointSymmetricPolytope`    | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `APolygon`                   | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `APolytope`                  |  ⊆  | (⊆) | [⊆] | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `ASingleton`                 |{⊆},∩=∅|(⊆,∩=∅)|**!⊆!**,**!∩=∅!**|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|[⊆,∩=∅]|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|((⊆,∩=∅))|(⊆,∩=∅)|((⊆,∩=∅))|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|
|                              |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| **Basic set types**          |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Ball1`                      | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Ball2`                      |     |     | (⊆) |     |      |     |     |{⊆},(∩=∅)|   |⊆,∩=∅|     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `BallInf`                    | (⊆) | (⊆) |(⊆,∩=∅)|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Ballp`                      |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Ellipsoid`                  |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `EmptySet`                   |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `HPolygon`/`HPolygonOpt`     | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `HPolytope`                  | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Hyperrectangle`             | (⊆) | (⊆) |(⊆,∩=∅)|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `Singleton`                  |(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|
| `VPolygon`                   | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
| `ZeroSet`                    |(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|(⊆,∩=∅)|
| `Zonotope`                   | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |
|                              |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| **Lazy set operation types** |     |     |     |     |      |     |     |       |     |     |     |     |     |     |     |     |      |       |     |       |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `CartesianProduct`           |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `CartesianProductArray`      |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ConvexHull`                 |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ConvexHullArray`            |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ExponentialMap`             |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `ExponentialProjectionMap`   |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `HalfSpace`                  |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Hyperplane`                 |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `Intersection`               |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `LinearMap`                  |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `MinkowskiSum`               |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `MinkowskiSumArray`          |     |     | (⊆) |     |      |     |     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |      |(⊆,∩=∅)|     |(⊆,∩=∅)|     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| `SymmetricIntervalHull`      | (⊆) | (⊆) |(⊆,∩=∅)|(⊆)| (⊆)  | (⊆) | (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆)  |(⊆,∩=∅)| (⊆) |(⊆,∩=∅)| (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) | (⊆) |


### `⊆`

This function checks whether a set is a subset of another set.
It can optionally produce a witness if the subset relation does not hold.
The operator can be used in infix notation (`X ⊆ S`).
Alias: `issubset`

### `is_intersection_empty`

This function checks whether the intersection of two sets is empty.
It can optionally produce a witness if the intersection is nonempty.
