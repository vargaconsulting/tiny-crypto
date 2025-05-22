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
ECPoint(curve::C) where {C<:Curve} = ECPoint(typeof(curve.G.point.x), curve)
Infinity(curve::C) where {C<:Curve} = ECPoint(typeof(curve.G.point.x), curve)

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

function curve_points(curve::C) where {C<:Curve}
    T = typeof(curve.G.point.x)
    π = Int(T.parameters[2])
    points = Point{T}[]

    for x_raw in 0:π-1
        x = T(x_raw)
        rhs = curve_equation_rhs(x, curve)

        for y_raw in 0:π-1
            y = T(y_raw)
            lhs = curve_equation_lhs(Point{T}(x, y), curve)
            lhs == rhs && push!(points, Point{T}(x, y))
        end
    end

    push!(points, infinity(T))
    return points
end

function is_point_on_curve(P::Point{F}, curve::C) where {F<:Fp, C<:Curve}
    is_infinity(P) && return true
    return curve_equation_lhs(P, curve) == curve_equation_rhs(P.x, curve)
end
is_point_on_curve(P::ECPoint{F,C}) where {F<:Fp, C<:Curve} = is_point_on_curve(P.point, P.curve)
is_point_on_curve(P::ECPoint{F,C}, curve::C) where {F<:Fp, C<:Curve} = is_point_on_curve(P.point, curve)

function is_generator(G::Point{F}, curve::C, n::Integer) where {F<:Fp, C<:Curve}
    is_point_on_curve(G, curve) || return false
    return scalar_mult(n, G, curve) == infinity(F)
end
is_generator(P::ECPoint{T,C}, n::Integer) where {T<:Fp, C<:Curve} = is_generator(P.point, P.curve, n)
is_generator(P::ECPoint{T,C}, curve::C, n::Integer) where {T<:Fp, C<:Curve} = is_generator(P.point, curve, n)

function is_singular(::C) where {C<:Curve}
    error("Singularity test not implemented for curve type $(C)")
end

function point_add(P::Point{T}, Q::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    error("point_add not implemented for curve type $(typeof(curve))")
end

function scalar_mult(k::Integer, P::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    k < 0 && return scalar_mult(-k, Point{T}(P.x, -P.y), curve)
    R = infinity(T)
    N = P
    while k > 0
        if k & 1 != 0
            R = point_add(R, N, curve)
        end
        N = point_add(N, N, curve)
        k >>= 1
    end
    return R
end

function point_order(P::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    if is_infinity(P) return 1 end
    Q, E = P, curve_points(curve)
    N = length(E)  # total number of curve points including ∞
    for i in 1:N
        if is_infinity(Q) return i end
        Q = point_add(Q, P, curve)
    end
    return nothing
end
point_order(P::ECPoint{T,C}) where {T<:Fp, C<:Curve} = point_order(P.point, P.curve)

