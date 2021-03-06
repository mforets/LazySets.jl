export SymmetricIntervalHull,
       an_element

"""
    SymmetricIntervalHull{N, S<:LazySet{N}} <: AbstractHyperrectangle{N}

Type that represents the symmetric interval hull of a compact convex set.

### Fields

- `X`     -- compact convex set
- `cache` -- partial storage of already computed bounds, organized as mapping
    from dimension to tuples `(bound, valid)`, where `valid` is a flag
    indicating if the `bound` entry has been computed

### Notes

The symmetric interval hull can be computed with ``2n`` support vector queries
of unit vectors, where ``n`` is the dimension of the wrapped set (i.e., two
queries per dimension).
When asking for the support vector for a direction ``d``, one needs ``2k`` such
queries, where ``k`` is the number of non-zero entries in ``d``.

However, if one asks for many support vectors in a loop, the number of
computations may exceed ``2n``.
To be most efficient in such cases, this type stores the intermediately computed
bounds in the `cache` field.

The set `X` must be compact.
"""
struct SymmetricIntervalHull{N, S<:LazySet{N}} <: AbstractHyperrectangle{N}
    X::S
    cache::Vector{N}

    # default constructor that initializes cache
    function SymmetricIntervalHull(X::S) where {N, S<:LazySet{N}}
        @assert isbounded(X) "the symmetric interval hull is only defined " *
                             "for bounded sets"
        cache = fill(-one(N), dim(X))
        return new{N, S}(X, cache)
    end
end

isoperationtype(::Type{<:SymmetricIntervalHull}) = true
isconvextype(::Type{<:SymmetricIntervalHull}) = true

"""
    SymmetricIntervalHull(∅::EmptySet)

The symmetric interval hull of an empty set.

### Input

- `∅` -- an empty set

### Output

The empty set because it is absorbing for the symmetric interval hull.
"""
SymmetricIntervalHull(∅::EmptySet) = ∅


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(sih::SymmetricIntervalHull, i::Int)

Return the box radius of the symmetric interval hull of a convex set in a given
dimension.

### Input

- `sih` -- symmetric interval hull of a convex set
- `i`   -- dimension of interest

### Output

The radius in the given dimension.
If it was computed before, this is just a look-up, otherwise it requires two
support vector computations.
"""
function radius_hyperrectangle(sih::SymmetricIntervalHull, i::Int)
    return get_radius!(sih, i)
end

"""
    radius_hyperrectangle(sih::SymmetricIntervalHull)

Return the box radius of the symmetric interval hull of a convex set in every
dimension.

### Input

- `sih` -- symmetric interval hull of a convex set

### Output

The box radius of the symmetric interval hull of a convex set.

### Notes

This function computes the symmetric interval hull explicitly.
"""
function radius_hyperrectangle(sih::SymmetricIntervalHull)
    n = dim(sih)
    for i in 1:n
        get_radius!(sih, i, n)
    end
    return sih.cache
end

"""
    center(sih::SymmetricIntervalHull{N}, i::Int) where {N}

Return the center along a given dimension of a symmetric interval hull of a convex set.

### Input

- `sih` -- symmetric interval hull of a convex set
- `i`   -- dimension of interest

### Output

The center along a given dimension of the symmetric interval hull of a convex set.
"""
@inline function center(sih::SymmetricIntervalHull{N}, i::Int) where {N}
    @boundscheck _check_bounds(sih, i)
    return zero(N)
end


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(sih::SymmetricIntervalHull{N}) where {N}

Return the center of a symmetric interval hull of a convex set.

### Input

- `sih` -- symmetric interval hull of a convex set

### Output

The origin.
"""
function center(sih::SymmetricIntervalHull{N}) where {N}
    return zeros(N, dim(sih))
end


# --- LazySet interface functions ---


"""
    dim(sih::SymmetricIntervalHull)

Return the dimension of a symmetric interval hull of a convex set.

### Input

- `sih` -- symmetric interval hull of a convex set

### Output

The ambient dimension of the symmetric interval hull of a convex set.
"""
function dim(sih::SymmetricIntervalHull)
    return dim(sih.X)
end

"""
    σ(d::AbstractVector, sih::SymmetricIntervalHull)

Return the support vector of a symmetric interval hull of a convex set in a
given direction.

### Input

- `d`   -- direction
- `sih` -- symmetric interval hull of a convex set

### Output

The support vector of the symmetric interval hull of a convex set in the given
direction.
If the direction has norm zero, the origin is returned.

### Algorithm

For each non-zero entry in `d` we need to either look up the bound (if it has
been computed before) or compute it, in which case we store it for future
queries.
One such computation just asks for the support vector of the underlying set for
both the positive and negative unit vector in the respective dimension.
"""
function σ(d::AbstractVector, sih::SymmetricIntervalHull)
    n = length(d)
    @assert n == dim(sih) "cannot compute the support vector of a " *
        "$(dim(sih))-dimensional set along a vector of length $n"
    N = promote_type(eltype(d), eltype(sih))
    svec = similar(d)
    for i in eachindex(d)
        if d[i] == zero(N)
            svec[i] = zero(N)
        else
            svec[i] = sign(d[i]) * get_radius!(sih, i, n)
        end
    end
    return svec
end


# --- SymmetricIntervalHull functions ---


"""
    get_radius!(sih::SymmetricIntervalHull{N},
                i::Int,
                n::Int=dim(sih)) where {N}

Compute the radius of the symmetric interval hull of a convex set in a given
dimension.

### Input

- `sih` -- symmetric interval hull of a convex set
- `i`   -- dimension in which the radius should be computed
- `n`   -- (optional, default: `dim(sih)`) set dimension

### Output

The radius of the symmetric interval hull of a convex set in a given dimension.

### Algorithm

We ask for the support vector of the underlying set for both the positive and
negative unit vector in the dimension `i`.
"""
function get_radius!(sih::SymmetricIntervalHull{N},
                     i::Int,
                     n::Int=dim(sih)) where {N}
    if sih.cache[i] == -one(N)
        right_bound = σ(sparsevec([i], [one(N)], n), sih.X)
        left_bound = σ(sparsevec([i], [-one(N)], n), sih.X)
        sih.cache[i] = max(right_bound[i], abs(left_bound[i]))
    end
    return sih.cache[i]
end

# Faster support function calculation for sev and SymmetricIntervalHull
function ρ(d::SingleEntryVector, H::SymmetricIntervalHull)
    @assert length(d) == dim(H) "a $(d.n)-dimensional vector is " *
                                "incompatible with a $(dim(H))-dimensional set"
    return abs(d.v) * radius_hyperrectangle(H, d.i)
end

# Faster support vector calculation for sev and SymmetricIntervalHull
function σ(d::SingleEntryVector, sih::SymmetricIntervalHull)
    N = promote_type(eltype(d), eltype(sih))
    @assert length(d) == dim(sih) "cannot compute the support vector of a " *
                                  "$(dim(sih))-dimensional set along a vector of length $(d.n)"
    if d.v < zero(N)
        return SingleEntryVector(d.i, d.n, -radius_hyperrectangle(sih, d.i))
    else
        return SingleEntryVector(d.i, d.n, radius_hyperrectangle(sih, d.i))
    end
end
