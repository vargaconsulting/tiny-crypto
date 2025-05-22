struct Point{F<:Fp}
    x::Union{F, Nothing}
    y::Union{F, Nothing}
end
# Converting constructor: allows Point{S} -> Point{T}
Point{F}(p::Point) where {F<:Fp} = Point{F}(
    convert(Union{F, Nothing}, p.x),
    convert(Union{F, Nothing}, p.y)
)
# Equality
Base.:(==)(A::Point{F}, B::Point{F}) where {F<:Fp} = A.x == B.x && A.y == B.y

# Infinity
infinity(::Type{F}) where {F<:Fp} = Point{F}(nothing, nothing)
is_infinity(P::Point) = P.x === nothing && P.y === nothing

# Show
function Base.show(io::IO, P::Point{F}) where {F<:Fp}
    if is_infinity(P)
        print(io, "(∞,∞)")
    else
        print(io, "(", P.x, ",", P.y, ")")
    end
end

abstract type Curve end
struct ECPoint{T,C<:Curve}
    point::Point{T}
    curve::C # TODO: measure performance and replace this with reference if it makes sense
end
function ECPoint(::Type{T}, curve::C) where {T<:Fp, C<:Curve}
    return ECPoint{T,C}(infinity(T), curve)
end
ECPoint(curve::C) where {C<:Curve} = ECPoint(typeof(curve.G.x), curve)
Infinity(curve::C) where {C<:Curve} = ECPoint(typeof(curve.G.x), curve)

Base.:(==)(A::ECPoint{T,C}, B::ECPoint{T,C}) where {T,C<:Curve} =
    A.point == B.point && A.curve == B.curve

Base.show(io::IO, P::ECPoint{T,C}) where {T,C} =
    print(io, "$(P.point) ∈ ", typeof(P.curve))

Base.:+(A::ECPoint{T,C}, B::ECPoint{T,C}) where {T,C<:Curve} = begin
    A.curve == B.curve || error("ECPoint curve mismatch")
    ECPoint(point_add(A.point, B.point, A.curve), A.curve)
end

Base.:*(k::Integer, A::ECPoint{T,C}) where {T,C<:Curve} = ECPoint(scalar_mult(k, A.point, A.curve), A.curve)
Base.:-(A::ECPoint{T,C}) where {T,C<:Curve} = is_infinity(A.point) ? A : ECPoint(Point{T}(A.point.x, -A.point.y), A.curve)
Base.:-(A::ECPoint{T,C}, B::ECPoint{T,C}) where {T,C<:Curve} = A + (-B)

is_infinity(P::ECPoint) = is_infinity(P.point)
