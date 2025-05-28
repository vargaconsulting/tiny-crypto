using ..Macroes
using ..Utils: is_prime
using ..Field: Fp, 𝔽ₚ

mutable struct Weierstrass{F<:Fp} <: Curve
    @define(Int, π, a, b, order, cofactor)
    G::ECPoint{F, Weierstrass{F}}

    function Weierstrass{F}(a::Integer, b::Integer, order::Integer, cofactor::Integer, Gxy::Tuple{Integer,Integer}) where {F<:Fp}
        π = Int(F.parameters[2])
        self = new{F}(π, a, b, order, cofactor) # Incomplete initialization step
        self.G = ECPoint(Point{F}(F(Gxy[1]), F(Gxy[2])), self) # Fully initialize G
        return self
    end
end

function Weierstrass(π::T, a::Integer, b::Integer, order::Integer, cofactor::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    F = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return Weierstrass{F}(a, b, order, cofactor, G)
end
curve_equation_rhs(P::Point{F}, curve::Weierstrass{F}) where {F<:Fp} = P.x^3 + curve.a * P.x + curve.b
curve_equation_lhs(P::Point{F}, curve::Weierstrass{F}) where {F<:Fp} = P.y^2
infinity(curve::Weierstrass{F}) where {F<:Fp} = Point{F}(nothing, nothing)
is_complete(curve::Weierstrass{F}) where {F<:Fp} = false

function is_singular(curve::Weierstrass{F}) where {F<:Fp}
    @attach(curve, a,b, π)
    return mod(4a^3 + 27b^2, π) == 0
end

function Base.show(io::IO, curve::Weierstrass{F}) where {F<:Fp}
    @attach(curve, π, a, b, order, cofactor, G)

    a_int, b_int, q, h = Int(a), Int(b), order, cofactor
    gx, gy = Int(G.point.x), Int(G.point.y)

    bold(x) = "\e[1m$(x)\e[22m"
    subscript(x) = join(Char(0x2080 + d) for d in reverse(digits(x)))

    terms = ["x³"]
    a_int != 0 && push!(terms, "$(bold(a_int))x")
    b_int != 0 && push!(terms, bold(b_int))
    eq_str = join(terms, " + ")

    field_str = "𝔽$(subscript(π))"
    curve_type = "Weierstrass{$field_str}"
    generator_str = "𝔾($(bold(gx)),$(bold(gy)))"

    print(io, "$curve_type: y² = $eq_str | $generator_str, q = $(bold(q)), h = $(bold(h)), #E = $(bold(q * h))")
end

function point_neg(P::Point{F}, curve::Weierstrass{F}) where {F<:Fp}
    @attach(P, x,y)
    is_infinity(P) && return P
    return Point{F}(x, -y)
end

function point_add(P::Point{F}, Q::Point{F}, curve::Weierstrass{F})::Point{F} where {F<:Fp}
    @attach(P, x,y)
    is_infinity(P) && return Q
    is_infinity(Q) && return P
    (P.x == Q.x && P.y != Q.y) && return infinity(curve)

    λ = if P != Q
        (Q.y - P.y) / (Q.x - P.x)
    else
        (F(3) * P.x^2 + curve.a) / (F(2) * P.y)
    end

    x3 = λ^2 - P.x - Q.x
    y3 = λ * (P.x - x3) - P.y
    return Point{F}(x3, y3)
end

## we will rework this later
function Weierstrass(Π::UnitRange, A::UnitRange, B::UnitRange, ::Val{T}=Val(UInt128)) where {T<:Unsigned}
    for π ∈ primes(Π)
        is_prime(π) || continue
        for a ∈ A, b ∈ B
            F = Fp{T, π}
            proposed = Weierstrass{F}(a, b, 0, 1, (0, 0))  # dummy order, dummy generator
            is_singular(proposed) && continue
            E = curve_points(proposed)
            N = length(E)
            
            is_prime(N) || continue
            
            for P in E
                if !is_infinity(P) && is_generator(P, proposed, N)
                    return Weierstrass{F}(a, b, N, 1, (Int(P.x), Int(P.y)))
                end
            end
        end
    end
    @info("No suitable curve found.")
end

