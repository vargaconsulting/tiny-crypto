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

# Show

abstract type Curve end
struct AffinePoint{T,C<:Curve}
    point::Point{T}
    curve::C
end
AffinePoint(curve::C) where {C<:Curve} = AffinePoint(typeof(curve.G.point.x), curve)

# Infinity
is_infinity(P::Point) = P.x === nothing && P.y === nothing
identity(curve::C) where {C<:Curve} = infinity(curve)

Base.zero(curve::C) where {C<:Curve} = AffinePoint(identity(curve), curve)
Base.identity(curve::C) where {C<:Curve} = AffinePoint(identity(curve), curve)
Base.:(==)(A::Point{F}, B::Point{F}) where {F<:Fp} = A.x == B.x && A.y == B.y
Base.:(==)(A::AffinePoint{T,C}, B::AffinePoint{T,C}) where {T,C<:Curve} = A.point == B.point && A.curve == B.curve
Base.:+(A::AffinePoint{T,C}, B::AffinePoint{T,C}) where {T,C<:Curve} = begin
    A.curve == B.curve || error("AffinePoint curve mismatch")
    AffinePoint(point_add(A.point, B.point, A.curve), A.curve)
end
Base.:*(k::Integer, A::AffinePoint{T,C}) where {T,C<:Curve} = is_infinity(A) || is_identity(A) ? A : AffinePoint(scalar_mult(k, A.point, A.curve), A.curve)

Base.:-(A::AffinePoint{T,C}, B::AffinePoint{T,C}) where {T,C<:Curve} = A + (-B)
point_neg(P::Point{T}, curve::C) where {T<:Fp, C<:Curve} = error("point_neg not implemented for curve $(typeof(curve))")
Base.:-(A::AffinePoint{T,C}) where {T<:Fp, C<:Curve} = is_infinity(A) ? A : AffinePoint(point_neg(A.point, A.curve), A.curve)

is_infinity(P::AffinePoint) = is_infinity(P.point)
is_identity(P::Point{F}, curve::C) where {F<:Fp, C<:Curve} = P == identity(curve)
is_identity(P::AffinePoint{F,C}) where {F<:Fp, C<:Curve} = is_identity(P.point, P.curve)
is_identity(P::AffinePoint) = is_identity(P.point, P.curve)
inverse(A::AffinePoint) = is_infinity(A) || is_identity(A) ? A : AffinePoint(point_neg(A.point, A.curve), A.curve)

function curve_points(curve::C) where {C<:Curve}
    T = typeof(curve.G.point.x)
    π = Int(T.parameters[2])
    points = Point{T}[]

    for x_raw in 0:π-1
        x = T(x_raw)
        for y_raw in 0:π-1
            y = T(y_raw)
            P = Point{T}(x, y)
            lhs = curve_equation_lhs(P, curve)
            rhs = curve_equation_rhs(P, curve)
            lhs == rhs && push!(points, P)
        end
    end
    identity_point = identity(curve)
    in(identity_point, points) || push!(points, identity_point)
    return points
end
function subgroup_points(curve::C) where {C<:Curve}
    G = curve.G
    points = Set([k * G for k in 1:curve.order - 1])
    return union(points, [AffinePoint(identity(curve), curve)])
end
function is_point_on_curve(P::Point{F}, curve::C) where {F<:Fp, C<:Curve}
    is_infinity(P) && return true
    return curve_equation_lhs(P, curve) == curve_equation_rhs(P, curve)
end
is_point_on_curve(P::AffinePoint{F,C}) where {F<:Fp, C<:Curve} = is_point_on_curve(P.point, P.curve)
is_point_on_curve(P::AffinePoint{F,C}, curve::C) where {F<:Fp, C<:Curve} = is_point_on_curve(P.point, curve)

function is_generator(G::Point{F}, curve::C, n::Integer) where {F<:Fp, C<:Curve}
    is_point_on_curve(G, curve) || return false
    return scalar_mult(n, G, curve) == identity(curve)
end
is_generator(P::AffinePoint{T,C}, n::Integer) where {T<:Fp, C<:Curve} = is_generator(P.point, P.curve, n)
is_generator(P::AffinePoint{T,C}, curve::C, n::Integer) where {T<:Fp, C<:Curve} = is_generator(P.point, curve, n)

function is_singular(::C) where {C<:Curve}
    error("Singularity test not implemented for curve type $(C)")
end

function point_add(P::Point{T}, Q::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    error("point_add not implemented for curve type $(typeof(curve))")
end
function point_add(P::AffinePoint{T,C}, Q::AffinePoint{T,C}) where {T<:Fp, C<:Curve}
    P.curve == Q.curve || error("Mismatched curves")
    return AffinePoint(point_add(P.point, Q.point, P.curve), P.curve)
end

function scalar_mult_(k::Integer, P::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    if k == 0 || is_infinity(P)
        return identity(curve)
    elseif k < 0
        return scalar_mult(-k, point_neg(P, curve), curve)
    end

    R = identity(curve)
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
function scalar_mult(k::Integer, P::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    if k == 0 || is_infinity(P)
        return identity(curve)
    elseif k < 0
        return scalar_mult(-k, point_neg(P, curve), curve)
    end

    R = identity(curve)
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

function has_order(P::Point{T},  q::Int, curve::C) where {T<:Fp, C<:Curve}
    is_infinity(P) && return false
    G = AffinePoint(P, curve)
    return is_identity(q * G)
end

function point_order(P::Point{T}, E::Vector{Point{T}}, curve::C) where {T<:Fp, C<:Curve}
    is_infinity(P) && return 1
    Q, max_order = P, length(E)
    for i in 1:max_order
        if is_identity(Q, curve)
            return i
        end
        Q = point_add(Q, P, curve)
    end
    return curve.order
end

function point_order(P::Point{T}, curve::C) where {T<:Fp, C<:Curve}
    (is_infinity(P) || is_identity(P, curve)) && return 1
    Q = P
    for n in 2:curve.order
        Q = point_add(Q, P, curve)
        is_identity(Q, curve) && return n
    end
    return 1
end
point_order(P::AffinePoint{T,C}) where {T<:Fp, C<:Curve} = point_order(P.point, P.curve)

