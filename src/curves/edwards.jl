using ..Macroes
using ..Utils: is_prime, is_square
using ..Field: Fp, 𝔽ₚ

mutable struct TwistedEdwards{F<:Fp} <: Curve
    @define(Int, π, a, d, order, cofactor)
    G::ECPoint{F, TwistedEdwards{F}}

    function TwistedEdwards{F}(a::Integer, d::Integer, order::Integer, cofactor::Integer, Gxy::Tuple{Integer,Integer}) where {F<:Fp}
        π = Int(F.parameters[2])
        self = new{F}(π, a, d, order, cofactor)
        self.G = ECPoint(Point{F}(F(Gxy[1]), F(Gxy[2])), self)
        return self
    end
end

function TwistedEdwards(π::T, a::Integer, d::Integer, order::Integer, cofactor::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    F = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return TwistedEdwards{F}(a, d, order, cofactor, G)
end
curve_equation_rhs(P::Point{F}, curve::TwistedEdwards{F}) where {F<:Fp} = F(1) + curve.d * P.x^2 * P.y^2
curve_equation_lhs(P::Point{F}, curve::TwistedEdwards{F}) where {F<:Fp} = curve.a * P.x^2 + P.y^2
infinity(curve::TwistedEdwards{F}) where {F<:Fp} = Point{F}(zero(F), one(F))

# some curves are incomplete, leading to 0 denominator in the point addition
function is_complete(curve::TwistedEdwards{F}) where {F<:Fp}
    @attach(curve, a, d, π)
    return is_square(a, π) && !is_square(d, π)
end

function is_singular(curve::TwistedEdwards{F}) where {F<:Fp}
    @attach(curve, a, d, π)
    return a == 0 || d == 0 || a == d || !is_complete(curve)
end

function Base.show(io::IO, curve::TwistedEdwards{F}) where {F<:Fp}
    @attach(curve, π, a, d, order, cofactor, G)

    a_int, d_int, q, h = Int(a), Int(d), order, cofactor
    gx, gy = Int(G.point.x), Int(G.point.y)

    bold(x) = "\e[1m$(x)\e[22m"
    subscript(x) = join(Char(0x2080 + d) for d in reverse(digits(x)))

    eq_str = "$(bold(a_int))x² + y² = 1 + $(bold(d_int))x²y²"
    field_str = "𝔽$(subscript(π))"
    curve_type = "TwistedEdwards{$field_str}"
    generator_str = "𝔾($(bold(gx)),$(bold(gy)))"

    print(io, "$curve_type: $eq_str | $generator_str, q = $(bold(q)), h = $(bold(h)), #E = $(bold(q * h))")
end

function point_neg(P::Point{F}, curve::TwistedEdwards{F}) where {F<:Fp}
    @attach(P, x,y)
    is_infinity(P) && return P           # for (nothing, nothing)
    is_identity(P, curve) && return P    # for (0, 1)
    return Point{F}(-x, y)
end

function point_add(P::Point{F}, Q::Point{F}, curve::TwistedEdwards{F})::Point{F} where {F<:Fp}
    (is_infinity(P) || is_identity(P, curve)) && return Q
    (is_infinity(Q) || is_identity(Q, curve)) && return P

    x₁, y₁, x₂, y₂, 𝟙𝔽, a,d = P.x, P.y, Q.x, Q.y, F(1),curve.a, curve.d
    
    x = (x₁ * y₂ + y₁ * x₂) / (𝟙𝔽 + d * x₁ * x₂ * y₁ * y₂)
    y = (y₁ * y₂ - a * x₁ * x₂) / (𝟙𝔽 - d * x₁ * x₂ * y₁ * y₂)

    return Point{F}(x, y)
end

function TwistedEdwards(Π::UnitRange, A::UnitRange, D::UnitRange, max_cofactor::Int = 8, ::Val{T} = Val(UInt128)) where {T<:Unsigned}
    best_curve, best_order = nothing, 0
    for π in primes(Π)
        for a in filter(<(π), A), d in filter(<(π), D)
            F = Fp{T, π}
            proposed = TwistedEdwards{F}(a, d, 0,1, (0, 0))
            is_singular(proposed) && continue

            E = curve_points(proposed)
            proposed.order = N = length(E)
            for h in 2:8
                q = div(N, h)
                (N % h == 0 && is_prime(q)) || continue
                for P ∈ E[2:end]
                    Q = scalar_mult(h, P, proposed)
                    is_identity(Q, proposed) && continue
                    order = point_order(Q, proposed)
                    #is_generator(P, proposed, order) || continue
                    if q > best_order
                        best_order = q
                        best_curve = TwistedEdwards{F}(a, d, q, h, (Int(Q.x), Int(Q.y)))
                        @info("New best $best_curve")
                    end
                end    
            end
        end
    end
    isnothing(best_curve) ? @info("No suitable Twisted Edwards curve with generator found.") :
                            @info("Returning best effort curve with subgroup order $best_order")
    return best_curve
end
