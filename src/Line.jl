export Line

"""
    Line{N<:Real} <: LazySet{N}

Type that represents a line in 2D of the form ``a⋅x = b`` (i.e., a special case
of a `Hyperplane`).

### Fields

- `a` -- normal direction
- `b` -- constraint

### Examples

The line ``y = -x + 1``:

```jldoctest
julia> Line([1., 1.], 1.)
LazySets.Line{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)
```
"""
struct Line{N<:Real, V<:AbstractVector{N}} <: LazySet{N}
    a::V
    b::N

    # default constructor with length constraint
    function Line{N, V}(a::V, b::N) where {N<:Real, V<:AbstractVector{N}}
        @assert length(a) == 2 "lines must be two-dimensional"
        return new{N, V}(a, b)
    end
end

# type-less convenience constructor
Line(a::V, b::N) where {N<:Real, V<:AbstractVector{N}} = Line{N, V}(a, b)

# constructor from a LinearConstraint
Line(c::LinearConstraint{N}) where {N<:Real} = Line(c.a, c.b)


# --- LazySet interface functions ---


"""
    dim(L::Line)::Int

Return the ambient dimension of a line.

### Input

- `L` -- line

### Output

The ambient dimension of the line, which is 2.
"""
function dim(L::Line)::Int
    return 2
end

"""
    σ(d::V, L::Line) where {N<:Real, V<:AbstractVector{N}}

Return the support vector of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

The support vector in the given direction, which is defined the same way as for
the more general `Hyperplane`.
"""
function σ(d::V, L::Line) where {N<:Real, V<:AbstractVector{N}}
    return σ(d, Hyperplane(L.a, L.b))
end
