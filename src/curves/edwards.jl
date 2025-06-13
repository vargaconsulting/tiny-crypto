using ..Macroes
using ..Utils: is_prime, is_square
using ..Field: Fp, 𝔽ₚ

mutable struct Edwards{F<:Fp} <: Curve
    @define(Int, π, a, d, order, cofactor)
    G::AffinePoint{F, Edwards{F}}

    function Edwards{F}(a::Integer, d::Integer, order::Integer, cofactor::Integer, Gxy::Tuple{Integer,Integer}) where {F<:Fp}
        π = Int(F.parameters[2])
        self = new{F}(π, a, d, order, cofactor)
        self.G = AffinePoint(Point{F}(F(Gxy[1]), F(Gxy[2])), self)
        return self
    end
end

function Edwards(π::T, d::Integer, order::Integer, cofactor::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    F = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return Edwards{F}(1, d, order, cofactor, G)
end
function TwistedEdwards(π::T, a::Integer, d::Integer, order::Integer, cofactor::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    F = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return Edwards{F}(a, d, order, cofactor, G)
end
curve_equation_rhs(P::Point{F}, curve::Edwards{F}) where {F<:Fp} = F(1) + curve.d * P.x^2 * P.y^2
curve_equation_lhs(P::Point{F}, curve::Edwards{F}) where {F<:Fp} = curve.a * P.x^2 + P.y^2
infinity(curve::Edwards{F}) where {F<:Fp} = Point{F}(zero(F), one(F))

# some curves are incomplete, leading to 0 denominator in the point addition
function is_complete(curve::Edwards{F}) where {F<:Fp}
    @attach(curve, a, d, π)
    return is_square(a, π) && !is_square(d, π)
end

function is_singular(curve::Edwards{F}) where {F<:Fp}
    @attach(curve, a, d, π)
    return a == 0 || d == 0 || a == d || !is_complete(curve)
end

function Base.show(io::IO, curve::Edwards{F}) where {F<:Fp}
    @attach(curve, π, a, d, order, cofactor, G)
    a_int, d_int = Int(a), Int(d)
    gx, gy = Int(G.point.x), Int(G.point.y)
    subscript = join(Char(0x2080 + d) for d in reverse(digits(π)))
    field_str = "𝔽$subscript"
    curve_name = a == 1 ? "Edwards" : "TwistedEdwards"
    
    bold(x) = "\e[1m$(x)\e[22m"
    eq_str = a == F(1) ? 
        "x² + y² = 1 + $(bold(d_int))x²y²" :
        "$(bold(a_int))x² + y² = 1 + $(bold(d_int))x²y²"

    println(io, "$(curve_name){$field_str}: $eq_str | 𝔾($(bold(gx)),$(bold(gy))), q = $(bold(order)), h = $(bold(cofactor)), #E = $(bold(order*cofactor))")
end

function point_neg(P::Point{F}, curve::Edwards{F}) where {F<:Fp}
    @attach(P, x,y)
    is_infinity(P) && return P           # for (nothing, nothing)
    is_identity(P, curve) && return P    # for (0, 1)
    return Point{F}(-x, y)
end

function point_add(P::Point{F}, Q::Point{F}, curve::Edwards{F})::Point{F} where {F<:Fp}
    (is_infinity(P) || is_identity(P, curve)) && return Q
    (is_infinity(Q) || is_identity(Q, curve)) && return P

    x₁, y₁, x₂, y₂, 𝟙𝔽, a,d = P.x, P.y, Q.x, Q.y, F(1),curve.a, curve.d
    
    x = (x₁ * y₂ + y₁ * x₂) / (𝟙𝔽 + d * x₁ * x₂ * y₁ * y₂)
    y = (y₁ * y₂ - a * x₁ * x₂) / (𝟙𝔽 - d * x₁ * x₂ * y₁ * y₂)

    return Point{F}(x, y)
end


function search_twisted_edwards(Π::AbstractVector{<:Integer}, A::AbstractVector{<:Integer}, D::AbstractVector{<:Integer}, 
    max_cofactor::Int = 8, ::Val{T} = Val(UInt128)) where {T<:Unsigned}
    best_curve, best_order = nothing, 0
    for π in primes(Π)
        for a in filter(<(π), A), d in filter(<(π), D)
            F = Fp{T, π}
            proposed = Edwards{F}(a, d, 0,1, (0, 0))
            is_singular(proposed) && continue

            E = curve_points(proposed)
            proposed.order = N = length(E)
            for h in 2:max_cofactor
                q = div(N, h)
                (N % h == 0 && is_prime(q)) || continue
                for P ∈ E[2:end]
                    Q = scalar_mult(h, P, proposed)
                    is_identity(Q, proposed) && continue
                    order = point_order(Q, proposed)
                    #is_generator(P, proposed, order) || continue
                    if q > best_order
                        best_order = q
                        best_curve = Edwards{F}(a, d, q, h, (Int(Q.x), Int(Q.y)))
                        @info("New best $best_curve")
                    end
                end    
            end
        end
    end
    isnothing(best_curve) ? @info("No suitable Edwards curve with generator found.") :
                            @info("Returning best effort curve with subgroup order $best_order")
    return best_curve
end

TwistedEdwards(Π, A, D, max_cofactor = 8, v::Val{T} = Val(UInt128)) where {T<:Unsigned} =
    search_twisted_edwards(Π, filter(≥(2), A), D, max_cofactor, v)

Edwards(Π, D; max_cofactor = 8, v::Val{T} = Val(UInt128)) where {T<:Unsigned} =
    search_twisted_edwards(Π, 1:1, D, max_cofactor, v)


