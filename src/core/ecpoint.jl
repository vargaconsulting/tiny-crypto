struct Point{T} # Parametric type for Point
    x::Union{T, Nothing}
    y::Union{T, Nothing}
end
# Converting constructor: allows Point{S} -> Point{T}
Point{T}(p::Point) where {T} = Point{T}(convert(Union{T, Nothing}, p.x), convert(Union{T, Nothing}, p.y))
Base.:(==)(A::Point{T}, B::Point{T}) where {T} = A.x == B.x && A.y == B.y

Base.convert(::Type{Point{T}}, p::Point) where {T} = Point{T}(
    convert(Union{T, Nothing}, p.x),
    convert(Union{T, Nothing}, p.y)
)

infinity(::Type{T}) where {T} = Point{T}(nothing, nothing)
is_infinity(P::Point) = P.x === nothing || P.y === nothing

const ∞ = Point{Int}(nothing, nothing)

function Base.show(io::IO, P::Point{T}) where T
    if is_infinity(P)
        print(io, "(∞,∞)")
    else
        print(io, "(", P.x, ",", P.y, ")")
    end
end

abstract type Curve end

struct ECPoint{T,C<:Curve}
    point::Point{T}
    # TODO: can I drop in replace this with the reference idea:  curve::Base.RefValue{C}  # or just a `const CurveRef = Ref{Curve}`
    curve::C 
end
function Base.:(==)(A::ECPoint{T,C}, B::ECPoint{T,C}) where {T,C<:Curve}
    return A.point == B.point && A.curve == B.curve
end

function Base.show(io::IO, P::ECPoint{T,C}) where {T,C}
    print(io, "$(P.point) ∈ ", typeof(P.curve))
end
function Base.:+(A::ECPoint{T,C}, B::ECPoint{T,C}) where {T,C<:Curve}
    A.curve == B.curve || error("ECPoint curve mismatch")
    R = point_add(A.point, B.point, A.curve)
    return ECPoint{T,C}(R, A.curve)
end

function Base.:*(k::Integer, A::ECPoint{T,C}) where {T,C<:Curve}
    R = scalar_mult(k, A.point, A.curve)
    return ECPoint{T,C}(R, A.curve)
end

function Base.:-(A::ECPoint{T,C}) where {T,C<:Curve}
    if is_infinity(A.point)
        return A
    end
    neg_point = Point{T}(A.point.x, (-A.point.y) % A.curve.π)
    return ECPoint{T,C}(neg_point, A.curve)
end


Base.:-(A::ECPoint{T,C}, B::ECPoint{T,C}) where {T,C<:Curve} = A + (-B)
