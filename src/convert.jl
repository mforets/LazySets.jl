import Base.convert

#= conversion between set types =#

# convert methods for identity (no-ops)
for T in subtypes(LazySet, true)
    @eval begin
        Base.convert(::Type{$T}, X::$T) = X
    end
end

"""
    convert(T::Type{HPOLYGON}, P::VPolygon) where {HPOLYGON<:AbstractHPolygon}

Converts a polygon in vertex representation to a polygon in constraint representation.

### Input

- `HPOLYGON` -- type used for dispatch
- `P`        -- polygon in vertex representation

### Output

A polygon in constraint representation.
"""
function convert(T::Type{HPOLYGON}, P::VPolygon) where {HPOLYGON<:AbstractHPolygon}
    return tohrep(P, T)
end

"""
    convert(::Type{HPOLYGON}, X::LazySet; [check_boundedness]::Bool=true) where {HPOLYGON<:AbstractHPolygon}

Converts a polyhedral set to a polygon in vertex representation.

### Input

- `PT` -- type used for dispatch
- `X`  -- set
- `check_boundedness` -- (optional, default `true`) if `true` check whether the
                         set `X` is bounded before creating the polygon

### Output

A polygon in constraint representation.

### Algorithm

We compute the list of constraints of `X`, then instantiate the polygon.
"""
function convert(::Type{HPOLYGON}, X::LazySet; check_boundedness::Bool=true) where {HPOLYGON<:AbstractHPolygon}
    @assert dim(X) == 2 "set must be two-dimensional for conversion, but it is of dimension $(dim(X))"
    PT = basetype(HPOLYGON)
    if check_boundedness && !isbounded(X)
        throw(ArgumentError("expected a bounded set for conversion to `HPolygon`"))
    end
    return PT(constraints_list(X))
end

"""
    convert(::Type{VPolygon}, P::AbstractHPolygon)

Converts a polygon in constraint representation to a polygon in vertex representation.

### Input

- `VPolygon` -- type used for dispatch
- `P`        -- polygon in constraint representation

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon}, P::AbstractHPolygon)
    return tovrep(P)
end

"""
    convert(::Type{HPolytope}, P::AbstractHPolygon)

Convert from polygon in H-representation to polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polygon

### Output

The polygon represented as 2D polytope.
"""
function convert(::Type{HPolytope}, P::AbstractHPolygon)
    return HPolytope(P.constraints)
end

"""
    convert(::Type{VPolytope}, P::AbstractPolytope)

Convert polytopic type to polytope in V-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The set `P` represented as a `VPolytope`.
"""
function convert(::Type{VPolytope}, P::AbstractPolytope)
    return VPolytope(vertices_list(P))
end

"""
    convert(::Type{VPolygon}, P::AbstractPolytope)

Convert polytopic set to polygon in V-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The 2D polytope represented as a polygon.
"""
function convert(::Type{VPolygon}, P::AbstractPolytope)
    @assert dim(P) == 2 "polytope must be two-dimensional for conversion"
    return VPolygon(vertices_list(P))
end

"""
    convert(::Type{VPolygon}, X::LazySet)

Generic conversion to polygon in vertex representation.

### Input

- `type` -- target type
- `X`    -- set

### Output

The 2D set represented as a polygon.

### Algorithm

We compute the list of vertices of `X` and wrap the result in a polygon in
vertex representation, which guarantees that the vertices are sorted in
counter-clockwise fashion.
"""
function convert(::Type{VPolygon}, X::LazySet)
    @assert dim(X) == 2 "set must be two-dimensional for conversion"
    return VPolygon(vertices_list(X))
end

"""
    convert(::Type{HPolytope}, P::AbstractPolytope)

Convert a polytopic set to a polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The given polytope represented as a polytope in constraint representation.

### Algorithm

First the list of constraints of `P` is computed, then the corresponding
`HPolytope` is created.
"""
function convert(::Type{HPolytope}, P::AbstractPolytope)
    return HPolytope(constraints_list(P))
end

"""
    convert(::Type{HPolyhedron}, P::AbstractPolytope)

Convert a polytopic set to a polyhedron in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The given polytope represented as a polyhedron in constraint representation.
"""
function convert(::Type{HPolyhedron}, P::AbstractPolytope)
    return HPolyhedron(constraints_list(P))
end

"""
    convert(::Type{HPolytope}, P::VPolytope)

Convert from polytope in V-representation to polytope in H-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The polytope in the dual representation.

### Algorithm

The `tohrep` function is invoked. It requires the `Polyhedra` package.
"""
function convert(::Type{HPolytope}, P::VPolytope)
    return tohrep(P)
end

"""
    convert(::Type{VPolytope}, X::LazySet; [prune]::Bool=true)

Generic conversion to polytope in vertex representation.

### Input

- `type`  -- target type
- `X`     -- set
- `prune` -- (optional, default: `true`) option to remove redundant vertices

### Output

The given set represented as a polytope in vertex representation.

### Algorithm

We compute the list of vertices of `X` and wrap the result in a polytope in
vertex representation, `VPolytope`. Use the option `prune` to select whether or not
to remove redundant vertices before constructing the polytope.
"""
function convert(::Type{VPolytope}, X::LazySet; prune::Bool=true)
    return VPolytope(vertices_list(X, prune=prune))
end

"""
    convert(::Type{VPolytope}, P::HPolytope)

Convert from polytope in H-representation to polytope in V-representation.

### Input

- `type` -- target type
- `P`    -- source polytope

### Output

The polytope in the dual representation.

### Algorithm

The `tovrep` function is invoked. It requires the `Polyhedra` package.
"""
function convert(::Type{VPolytope}, P::HPolytope)
    return tovrep(P)
end

"""
    convert(::Type{HPOLYGON}, P::HPolytope{N, VN};
            prune::Bool=true) where {N<:Real, VN<:AbstractVector{N}, HPOLYGON<:AbstractHPolygon}

Convert from 2D polytope in H-representation to polygon in H-representation.

### Input

- `type`  -- target type
- `P`     -- source polytope (must be 2D)
- `prune` -- (optional, default: `true`) flag for removing redundant
             constraints in the end
### Output

The 2D polytope represented as polygon.
"""
function convert(::Type{HPOLYGON}, P::HPolytope{N, VN};
                 prune::Bool=true) where {N<:Real, VN<:AbstractVector{N}, HPOLYGON<:AbstractHPolygon}
    @assert dim(P) == 2 "polytope must be two-dimensional for conversion"
    H = HPOLYGON(Vector{LinearConstraint{N, VN}}())
    for ci in constraints_list(P)
        addconstraint!(H, ci; prune=prune)
    end
    return H
end

# fast conversion from a 2D hyperrectangular set to a zonotope
function _convert_2D(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
    c = center(H)
    rx = radius_hyperrectangle(H, 1)
    ry = radius_hyperrectangle(H, 2)
    G = _genmat_2D(c, rx, ry)
    return Zonotope(c, G)
end

@inline function _genmat_2D(c::AbstractVector{N}, rx, ry) where {N}
    flat_x = isapproxzero(rx)
    flat_y = isapproxzero(ry)
    ncols = !flat_x + !flat_y
    G = Matrix{N}(undef, 2, ncols)
    if !flat_x
        @inbounds begin G[1] = rx; G[2] = zero(N) end
        if !flat_y
            @inbounds begin G[3] = zero(N); G[4] = ry end
        end
    elseif !flat_y
        @inbounds begin G[1] = zero(N); G[2] = ry end
    end
    return G
end

function load_genmat_2D_static()
return quote
    @inline function _genmat_2D(c::SVector{L, N}, rx, ry) where {L, N}
        flat_x = isapproxzero(rx)
        flat_y = isapproxzero(ry)
        if !flat_x && !flat_y
            G = SMatrix{2, 2, N, 4}(rx, zero(N), zero(N), ry)
        elseif !flat_x && flat_y
            G = SMatrix{2, 1, N, 2}(rx, zero(N))
        elseif flat_x && !flat_y
            G = SMatrix{2, 1, N, 2}(zero(N), ry)
        else
            G = SMatrix{2, 0, N, 0}()
        end
        return G
    end

    # this function is type-stable but doesn't prune the generators according
    # to flat dimensions of H
    function _convert_2D_static(::Type{Zonotope}, H::AbstractHyperrectangle{N}) where {N}
        c = center(H)
        rx = radius_hyperrectangle(H, 1)
        ry = radius_hyperrectangle(H, 2)
        G = SMatrix{2, 2, N, 4}(rx, zero(N), zero(N), ry)
        return Zonotope(c, G)
    end

    function _convert_static(::Type{Zonotope}, H::Hyperrectangle{N, <:SVector, <:SVector}) where {N}
        return Zonotope(center(H), _genmat_static(H))
    end
end end  # quote / load_genmat_2D_static

function convert(::Type{Zonotope}, S::Singleton{N, VN}) where {N, VN<:AbstractVector{N}}
    MT = LazySets.Arrays._matrix_type(VN)
    zero_genmat = MT(undef, dim(S), 0)
    return Zonotope(element(S), zero_genmat)
end

"""
    convert(::Type{Zonotope}, Z::AbstractZonotope)

Converts a zonotopic set to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `H`        -- zonotopic set

### Output

A zonotope.
"""
function convert(::Type{Zonotope}, Z::AbstractZonotope)
    return _convert_zonotope_fallback(Z)
end

function convert(::Type{Zonotope}, H::AbstractHyperrectangle)
    dim(H) == 2 && return _convert_2D(Zonotope, H)
    return _convert_zonotope_fallback(H)
end

_convert_zonotope_fallback(Z) = Zonotope(center(Z), genmat(Z))

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}) where {N<:Real,
        HN1<:AbstractHyperrectangle{N}, HN2<:AbstractHyperrectangle{N}}

Converts the cartesian product of two hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `cp`       -- cartesian product of two hyperrectangular sets

### Output

This method falls back to the conversion of the cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, HN1, HN2}
                ) where {N<:Real, HN1<:AbstractHyperrectangle{N},
                         HN2<:AbstractHyperrectangle{N}}
    return convert(Zonotope, convert(Hyperrectangle, cp))
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN})
        where {N<:Real, HN<:AbstractHyperrectangle{N}}

Converts the cartesian product array of hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `cpa`      -- cartesian product array of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method falls back to the conversion of the cartesian product to a single
hyperrectangle, and then from a hyperrectangle to a zonotope.
"""
function convert(::Type{Zonotope}, cpa::CartesianProductArray{N, HN}
                ) where {N<:Real, HN<:AbstractHyperrectangle{N}}
    return convert(Zonotope, convert(Hyperrectangle, cpa))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, ZN}
           ) where {N, ZN<:AbstractZonotope{N}}

Converts the lazy linear map of a zonotopic set to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `S`        -- linear map of a zonotopic set

### Output

A zonotope.

### Algorithm

This method first applies the (concrete) linear map to the zonotopic set and
then converts the result to a `Zonotope` type.
"""
function convert(::Type{Zonotope}, S::LinearMap{N, ZN}
                ) where {N, ZN<:AbstractZonotope{N}}
    return convert(Zonotope, linear_map(S.M, S.X))
end

"""
    convert(::Type{Zonotope}, S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
           ) where {N, HN1<:AbstractHyperrectangle{N},
                    HN2<:AbstractHyperrectangle{N}}

Converts the lazy linear map of the cartesian product of two hyperrectangular
sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `S`        -- linear map of the cartesian product of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the cartesian product to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope},
                 S::LinearMap{N, CartesianProduct{N, HN1, HN2}}
                ) where {N, HN1<:AbstractHyperrectangle{N},
                         HN2<:AbstractHyperrectangle{N}}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{Zonotope},S::LinearMap{N, CartesianProductArray{N, HN}}
           ) where {N, HN<:AbstractHyperrectangle{N}}

Converts the lazy linear map of the cartesian product of a finite number of
hyperrectangular sets to a zonotope.

### Input

- `Zonotope` -- type, used for dispatch
- `S`        -- linear map of a `CartesianProductArray` of hyperrectangular sets

### Output

A zonotope.

### Algorithm

This method first converts the cartesian product array to a zonotope, and then
applies the (concrete) linear map to the zonotope.
"""
function convert(::Type{Zonotope}, S::LinearMap{N, CartesianProductArray{N, HN}}
                ) where {N, HN<:AbstractHyperrectangle{N}}
    return linear_map(S.M, convert(Zonotope, S.X))
end

"""
    convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
           ) where {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from singleton to polygon in H-representation.

### Input

- `type` -- target type
- `S`    -- singleton

### Output

A polygon in constraint representation with the minimal number of constraints
(three).
"""
function convert(::Type{HPOLYGON}, S::AbstractSingleton{N}
                ) where {N<:Real, HPOLYGON<:AbstractHPolygon}
    constraints_list = Vector{LinearConstraint{N, Vector{N}}}(undef, 3)
    o = one(N)
    z = zero(N)
    v = element(S)
    constraints_list[1] = LinearConstraint([o, o], v[1] + v[2])
    constraints_list[2] = LinearConstraint([-o, z], -v[1])
    constraints_list[3] = LinearConstraint([z, -o], -v[2])
    return HPOLYGON(constraints_list)
end

"""
    convert(::Type{HPOLYGON}, L::LineSegment{N}
          ) where {N<:Real, HPOLYGON<:AbstractHPolygon}

Convert from line segment to polygon in H-representation.

### Input

- `type`  -- target type
- `L`     -- line segment
- `prune` -- (optional, default: `false`) flag for removing redundant
             constraints in the end
### Output

A flat polygon in constraint representation with the minimal number of
constraints (four).
"""
function convert(::Type{HPOLYGON}, L::LineSegment{N}) where {N<:Real, HPOLYGON<:AbstractHPolygon}
    H = HPOLYGON{N}()
    c = halfspace_left(L.p, L.q)
    addconstraint!(H, c; prune=false)
    addconstraint!(H, LinearConstraint(-c.a, -c.b); prune=false)
    line_dir = L.q - L.p
    c = LinearConstraint(line_dir, dot(L.q, line_dir))
    addconstraint!(H, c; prune=false)
    line_dir = -line_dir
    addconstraint!(H, LinearConstraint(line_dir, dot(L.p, line_dir)); prune=false)
    return H
end

"""
    convert(::Type{Hyperrectangle}, H::AbstractHyperrectangle)

Convert a hyperrectangular set to a hyperrectangle.

### Input

- `Hyperrectangle` -- hyperrectangle type, used for dispatch
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

"""
    convert(::Type{Interval}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to an interval.

### Input

- `Interval` -- interval type, used for dispatch
- `H`        -- hyperrectangular set

### Output

An interval.

### Examples

```jldoctest
julia> convert(Interval, Hyperrectangle([0.5], [0.5]))
Interval{Float64,IntervalArithmetic.Interval{Float64}}([0, 1])
```
"""
function convert(::Type{Interval}, H::AbstractHyperrectangle)
    @assert dim(H) == 1 "can only convert a one-dimensional $(typeof(H)) to `Interval`"
    return Interval(low(H, 1), high(H, 1))
end

"""
    convert(::Type{Interval}, S::LazySet{N}) where {N<:Real}

Converts a convex set to an interval.

### Input

- `Interval` -- interval type, used for dispatch
- `S`        -- convex set

### Output

An interval.
"""
function convert(::Type{Interval}, S::LazySet{N}) where {N<:Real}
    @assert dim(S) == 1 "cannot convert a $(dim(S))-dimensional $(typeof(S)) to `Interval`"
    return Interval(-ρ(N[-1], S), ρ(N[1], S))
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, HN}) where {N<:Real, HN<:AbstractHyperrectangle{N}}

Converts the cartesian product of a finite number of hyperrectangular sets to
a single hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `S`              -- cartesian product array of hyperrectangular set

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N, HN}) where {N<:Real, HN<:AbstractHyperrectangle{N}}
     n = dim(cpa)
     c = Vector{N}(undef, n)
     r = Vector{N}(undef, n)
     i = 1
     @inbounds for block_set in array(cpa)
         j = i + dim(block_set) - 1
         c[i:j] = center(block_set)
         r[i:j] = radius_hyperrectangle(block_set)
         i = j + 1
     end
     return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle},
            cp::CartesianProduct{N, HN1, HN2}) where {N<:Real, HN1<:AbstractHyperrectangle{N}, HN2<:AbstractHyperrectangle{N}}

Converts the cartesian product of two hyperrectangular sets to a single hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `S`              -- cartesian product of two hyperrectangular sets

### Output

A hyperrectangle.

### Algorithm

The result is obtained by concatenating the center and radius of each hyperrectangle.
This implementation uses the `center` and `radius_hyperrectangle` methods of
`AbstractHyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 cp::CartesianProduct{N, HN1, HN2}) where {N<:Real, HN1<:AbstractHyperrectangle{N}, HN2<:AbstractHyperrectangle{N}}
     X, Y = cp.X, cp.Y
     c = vcat(center(X), center(Y))
     r = vcat(radius_hyperrectangle(X), radius_hyperrectangle(Y))
     return Hyperrectangle(c, r)
end

"""
    convert(::Type{Hyperrectangle},
            cpa::CartesianProductArray{N, Interval{N}}) where {N<:Real}

Converts the cartesian product of a finite number of intervals to a single
hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `S`              -- cartesian product array of intervals

### Output

A hyperrectangle.

### Algorithm

This implementation uses the `min` and `max` methods of `Interval` to reduce
the allocatons and improve performance (see LazySets#1143).
"""
function convert(::Type{Hyperrectangle},
                 cpa::CartesianProductArray{N, Interval{N}}) where {N<:Real}
     # since the sets are intervals, the dimension of cpa is its length
     n = length(array(cpa))
     l = Vector{N}(undef, n)
     h = Vector{N}(undef, n)
     @inbounds for (i, Ii) in enumerate(array(cpa))
         l[i] = min(Ii)
         h[i] = max(Ii)
     end
     return Hyperrectangle(low=l, high=h)
end

"""
    convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N<:Real}

Converts a two-dimensional hyperrectangle to the cartesian product of two
intervals.

### Input

- `CartesianProduct` -- type used for dispatch
- `H`                -- hyperrectangle

### Output

The cartesian product of two intervals.
"""
function convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N<:Real}
    @assert dim(H) == 2 "the hyperrectangle must be two-dimensional to convert it to " *
            "the cartesian product of two intervals, but it is $(dim(H))-dimensional; " *
            "consider converting it to a 'CartesianProductArray{$N, Interval{$N}}' instead"
    Ix = Interval(low(H, 1), high(H, 1))
    Iy =  Interval(low(H, 2), high(H, 2))
    return CartesianProduct(Ix, Iy)
end

"""
    convert(::Type{CartesianProductArray{N, Interval{N}}},
            H::AbstractHyperrectangle{N}) where {N<:Real}

Converts a hyperrectangle to the cartesian product array of intervals.

### Input

- `CartesianProductArray` -- type used for dispatch
- `H`                     -- hyperrectangle

### Output

The cartesian product of a finite number of intervals.
"""
function convert(::Type{CartesianProductArray{N, Interval{N}}},
                 H::AbstractHyperrectangle{N}) where {N<:Real}
    Iarray = [Interval(low(H, i), high(H, i)) for i in 1:dim(H)]
    return CartesianProductArray(Iarray)
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, ZN1, ZN2}
           ) where {N<:Real, ZN1<:AbstractZonotope{N}, ZN2<:AbstractZonotope{N}}

Converts the cartesian product of two zonotopes to a new zonotope.

### Input

- `Zonotope` -- type used for dispatch
- `S`        -- cartesian product of two zonotopes

### Output

A zonotope.

### Algorithm

The cartesian product is obtained by:

- Concatenating the centers of each input zonotope.
- Arranging the generators in block-diagional fashion, and filled with zeros
  in the off-diagonal; for this reason, the generator matrix of the returned
  zonotope is built as a sparse matrix.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, ZN1, ZN2}
                ) where {N<:Real, ZN1<:AbstractZonotope{N}, ZN2<:AbstractZonotope{N}}
    Z1, Z2 = cp.X, cp.Y
    c = vcat(center(Z1), center(Z2))
    G = blockdiag(sparse(genmat(Z1)), sparse(genmat(Z2)))
    return Zonotope(c, G)
end

"""
    convert(::Type{HPolytope}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to a polytope in constraint representation.

### Input

- `HPolytope` -- type used for dispatch
- `H`         -- hyperrectangular set

### Output

A polytope in constraint representation.
"""
function convert(::Type{HPolytope},
                 H::AbstractHyperrectangle)
    return HPolytope(constraints_list(H))
end

"""
    convert(::Type{HPOLYGON}, H::AbstractHyperrectangle) where
        {HPOLYGON<:AbstractHPolygon}

Converts a hyperrectangular set to a polygon in constraint representation.

### Input

- `HPOLYGON`  -- type used for dispatch
- `H`         -- hyperrectangular set

### Output

A polygon in constraint representation.
"""
function convert(X::Type{HPOLYGON},
                 H::AbstractHyperrectangle) where {HPOLYGON<:AbstractHPolygon}
    @assert dim(H) == 2 "cannot convert a $(dim(H))-dimensional " *
        "hyperrectangle into a two-dimensional polygon"
    return HPOLYGON(constraints_list(H))
end

"""
    convert(::Type{IntervalArithmetic.Interval}, x::Interval)

Converts a `LazySets` interval to an `Interval` from `IntervalArithmetic`.

### Input

- `Interval` -- type used for dispatch, from `IntervalArithmetic`
- `x`        -- interval (`LazySets.Interval`)

### Output

An `IntervalArithmetic.Interval`.
"""
function convert(::Type{IntervalArithmetic.Interval}, x::Interval)
    return IntervalArithmetic.interval(min(x), max(x))
end

"""
    convert(::Type{Interval}, x::IntervalArithmetic.Interval)

Converts an `Interval` from `IntervalArithmetic` to an interval in `LazySets`.

### Input

- `Interval` -- type used for dispatch
- `x`        -- interval (`IntervalArithmetic.Interval`)

### Output

A `LazySets.Interval`.
"""
function convert(::Type{Interval}, x::IntervalArithmetic.Interval)
    return Interval(IntervalArithmetic.inf(x), IntervalArithmetic.sup(x))
end

"""
    convert(::Type{IntervalArithmetic.IntervalBox}, H::AbstractHyperrectangle)

Converts a hyperrectangular set to an `IntervalBox` from `IntervalArithmetic`.

### Input

- `IntervalBox` -- type used for dispatch
- `H`           -- hyperrectangular set

### Output

An `IntervalBox`.
"""
function convert(::Type{IntervalArithmetic.IntervalBox}, H::AbstractHyperrectangle)
    return IntervalArithmetic.IntervalBox(IntervalArithmetic.interval.(low(H), high(H)))
end

"""
    convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)

Converts an `IntervalBox` from `IntervalArithmetic` to a hyperrectangular set.

### Input

- `Hyperrectangle` -- type used for dispatch
- `IB`             -- interval box

### Output

A `Hyperrectangle`.

### Notes

`IntervalArithmetic.IntervalBox` uses *static* vectors to store each component
interval, hence the resulting `Hyperrectangle` has its center and radius represented
as a static vector (`SArray`).
"""
function convert(::Type{Hyperrectangle}, IB::IntervalArithmetic.IntervalBox)
    low_IB = IntervalArithmetic.inf.(IB)
    high_IB = IntervalArithmetic.sup.(IB)
    return Hyperrectangle(low=low_IB, high=high_IB)
end

"""
    convert(::Type{Hyperrectangle}, r::Rectification{N, AH})
        where {N<:Real, AH<:AbstractHyperrectangle{N}}

Converts a rectification of a hyperrectangle to a hyperrectangle.

### Input

- `Hyperrectangle` -- type used for dispatch
- `r`              -- rectification of a hyperrectangle

### Output

A `Hyperrectangle`.
"""
function convert(::Type{Hyperrectangle},
                 r::Rectification{N, AH}) where {N<:Real,
                                                 AH<:AbstractHyperrectangle{N}}
    return rectify(r.X)
end

"""
    convert(::Type{Interval},
            r::Rectification{N, IN}) where {N<:Real, IN<:Interval{N}}

Converts a rectification of an interval to an interval.

### Input

- `Interval` -- type used for dispatch
- `r`        -- rectification of an interval

### Output

An `Interval`.
"""
function convert(::Type{Interval},
                 r::Rectification{N, IN}) where {N<:Real, IN<:Interval{N}}
    return Interval(rectify([min(r.X), max(r.X)]))
end

"""
    convert(::Type{VPolytope},
            X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}

Converts the convex hull array of singletons to a polytope in V-representation.

### Input

- `VPolytope` -- type used for dispatch
- `X`         -- convex hull array of singletons

### Output

A polytope in vertex representation.
"""
function convert(::Type{VPolytope},
                 X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
    return VPolytope(vertices_list(X))
end

"""
    convert(::Type{VPolygon},
            X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}

Converts the convex hull array of singletons to a polygon in V-representation.

### Input

- `VPolygon`  -- type used for dispatch
- `X`         -- convex hull array of singletons

### Output

A polygon in vertex representation.
"""
function convert(::Type{VPolygon},
                 X::ConvexHullArray{N, Singleton{N, VT}}) where {N, VT}
    @assert dim(X) == 2
    return VPolygon(vertices_list(X))
end

"""
    convert(::Type{MinkowskiSumArray},
            X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}

Converts the Minkowski sum of a Minkowski sum array to a Minkowski sum array.

### Input

- `MinkowskiSumArray`  -- type used for dispatch
- `X`                  -- Minkowski sum of a Minkowski sum array

### Output

A Minkowski sum array.
"""
function convert(::Type{MinkowskiSumArray},
                 X::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}
    return MinkowskiSumArray(vcat(X.X, X.Y.array))
end

"""
    convert(::Type{Interval}, x::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval{N}}

Converts the Minkowski sum of two intervals into an interval.

### Input

- `Interval` -- type used for dispatch
- `x`        -- Minkowski sum of a pair of intervals

### Output

An interval.
"""
function convert(::Type{Interval}, x::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval{N}}
    return minkowski_sum(x.X, x.Y)
end

function convert(::Type{HalfSpace{N, Vector{N}}}, hs::HalfSpace{N, <:AbstractVector{N}}) where {N}
    return HalfSpace(Vector(hs.a), hs.b)
end

function convert(::Type{Zonotope},
                 am::AbstractAffineMap{N, <:AbstractZonotope{N}}) where {N<:Real}
    Z1 = convert(Zonotope, linear_map(matrix(am), set(am)))
    translate!(Z1, vector(am))
    return Z1
end

"""
    convert(::Type{HParallelotope}, Z::AbstractZonotope{N}) where {N}

Converts a zonotopic set of order one into a parallelotope in constraint
representation.

### Input

- `HParallelotope` -- type used for dispatch
- `Z`              -- zonotopic set

### Output

A parallelotope in constraint representation.

### Notes

This function requires that the list of constraints of `Z` are obtained in
the particular order returned from the constraints list function of a `Zonotope`.
"""
function convert(::Type{HParallelotope}, Z::AbstractZonotope{N}) where {N}
    @assert order(Z) == 1 "cannot convert a zonotope that is not of order 1 to"*
                          " a parallelotope"
    n = dim(Z)

    Z = convert(Zonotope, Z)
    constraints = constraints_list(Z) # TODO refactor + use _constraints_list_zonotope

    D = Matrix{N}(undef, n, n)
    c = Vector{N}(undef, 2n)
    j = 1
    @inbounds for i in 1:n
        D[i, :] = constraints[j].a
        c[i] = constraints[j].b
        c[i+n] = constraints[j+1].b
        j += 2
    end
    return HParallelotope(D, c)
end

"""
    convert(::Type{Zonotope}, cp::CartesianProduct{N, AZ1, AZ2}) where
         {N, AZ1<:AbstractZonotope{N}, AZ2<:AbstractZonotope{N}}

Converts a cartesian product of two zonotopes into a zonotope.

### Input

- `Zonotope` -- type used for dispatch
- `cp`       -- cartesian product of two zonotopes

### Output

A zonotope.

### Notes

This implementation creates a `Zonotope` with sparse vector and matrix representation.
"""
function convert(::Type{Zonotope}, cp::CartesianProduct{N, AZ1, AZ2}) where
         {N, AZ1<:AbstractZonotope{N}, AZ2<:AbstractZonotope{N}}
    Z1, Z2 = cp.X, cp.Y
    c = vcat(center(Z1), center(Z2))
    G = blockdiag(sparse(genmat(Z1)), sparse(genmat(Z2)))
    return Zonotope(c, G)
end

"""
    convert(::Type{Zonotope}, cpa::CartesianProductArray{N, AZ})
        where {N, AZ<:AbstractZonotope{N}}

Converts a cartesian product array of zonotopes into a zonotope.

### Input

- `Zonotope` -- type used for dispatch
- `cpa`       -- cartesian product array of zonotopes

### Output

A zonotope.

### Notes

This function requires the use of `SparseArrays`.
"""
function convert(::Type{Zonotope}, cpa::CartesianProductArray{N, AZ}) where
        {N, AZ<:AbstractZonotope{N}}
    arr = array(cpa)
    c = reduce(vcat, center.(arr))
    G = reduce(blockdiag, sparse.(genmat.(arr)))
    return Zonotope(c, G)
end

convert(::Type{HPolytope}, P::HPolyhedron) = HPolytope(copy(constraints_list(P)))
convert(::Type{HPolyhedron}, P::HPolytope) = HPolyhedron(copy(constraints_list(P)))

for T in [HPolygon, HPolygonOpt, HPolytope, HPolyhedron]
@eval begin
    function convert(::Type{$T}, P::Intersection{N, S1, S2}) where {N, S1<:AbstractPolyhedron{N}, S2<:AbstractPolyhedron{N}}
        clist = vcat(constraints_list(P.X), constraints_list(P.Y))
        return ($T)(clist)
    end

    function convert(::Type{$T}, P::IntersectionArray{N, S}) where {N, S<:AbstractPolyhedron}
        clist = reduce(vcat, constraints_list.(array(P)))
        return ($T)(clist)
    end
end
end
