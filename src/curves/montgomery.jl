using ..Macroes
using ..Utils: is_prime, prime_factors
using ..Field: Fp, 𝔽ₚ

mutable struct Montgomery{F<:Fp} <: Curve
    @define(Int, π, B, A, order, cofactor)
    G::AffinePoint{F, Montgomery{F}}

    function Montgomery{F}(B::Integer, A::Integer, order::Integer, cofactor::Integer, Gxy::Tuple{Integer,Integer}) where {F<:Fp}
        π = Int(F.parameters[2])
        self = new{F}(π, B, A, order, cofactor)
        self.G = AffinePoint(Point{F}(F(Gxy[1]), F(Gxy[2])), self)
        return self
    end
end

function Montgomery(π::T, B::Integer, A::Integer, order::Integer, cofactor::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    F = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return Montgomery{F}(B, A, order, cofactor, G)
end

curve_equation_rhs(P::Point{F}, curve::Montgomery{F}) where {F<:Fp} = P.x^3 + curve.A * P.x^2 + P.x
curve_equation_lhs(P::Point{F}, curve::Montgomery{F}) where {F<:Fp} = curve.B * P.y^2
infinity(curve::Montgomery{F}) where {F<:Fp} = Point{F}(nothing, nothing)

function is_singular(curve::Montgomery{F}) where {F<:Fp}
    @attach(curve, A, B, π)
    return B == 0 || mod(A^2 - 4, π) == 0
end

function Base.show(io::IO, curve::Montgomery{F}) where {F<:Fp}
    @attach(curve, π, A, B, order, cofactor, G)
    A_int, B_int = Int(A), Int(B)
    gx, gy = Int(G.point.x), Int(G.point.y)
    subscript = join(Char(0x2080 + d) for d in reverse(digits(π)))
    field_str = "𝔽$subscript"

    bold(x) = "\e[1m$(x)\e[22m"
    lhs_str = "$(bold(B_int))y²"
    rhs_terms = ["x³"]
    A_int != 0 && push!(rhs_terms, "$(bold(A_int))x²")
    push!(rhs_terms, "x")
    rhs_str = join(rhs_terms, " + ")

    println(io, "Montgomery{$field_str}: $lhs_str = $rhs_str | 𝔾($(bold(gx)),$(bold(gy))), q = $(bold(order)), h = $(bold(cofactor)), #E = $(bold(order * cofactor))")
end

function point_neg(P::Point{F}, curve::Montgomery{F}) where {F<:Fp}
    @attach(P, x,y)
    is_infinity(P) && return P
    return Point{F}(x, -y)
end

function point_add(P::Point{F}, Q::Point{F}, curve::Montgomery{F})::Point{F} where {F<:Fp}
    (is_infinity(P) || is_identity(P, curve)) && return Q
    (is_infinity(Q) || is_identity(Q, curve)) && return P

    x₁, y₁, x₂, y₂,  A, B = P.x, P.y, Q.x, Q.y, curve.A, curve.B
    𝟙𝔽, 𝟚𝔽, 𝟛𝔽 = F(1), F(2), F(3)
    if P == Q
        iszero(y₁) && return infinity(curve)
        λ = (𝟛𝔽 * x₁^2 + 𝟚𝔽 * A * x₁ + 𝟙𝔽) / (𝟚𝔽 * B * y₁)
    elseif x₁ == x₂
        return infinity(curve)
    else
        λ = (y₂ - y₁) / (x₂ - x₁)
    end

    x = B * λ^2 - A - x₁ - x₂
    y = λ * (x₁ - x) - y₁
    return Point{F}(x, y)
end


function Montgomery(Π::AbstractVector{<:Integer}, B::AbstractVector{<:Integer}, A::AbstractVector{<:Integer},
    max_cofactor::Int = 8, ::Val{T} = Val(UInt128)) where {T<:Unsigned}
    
    best_curve, best_order = nothing, 0
    for π in primes(Π)
        for b in filter(<(π), B), a in filter(<(π), A)
            F = Fp{T, π}
            proposed = Montgomery{F}(b, a, 0, 1, (0, 0))
            is_singular(proposed) && continue

            E = curve_points(proposed)
            proposed.order = N = length(E)

            for h in 2:max_cofactor
                q = div(N, h)
                (N % h == 0 && is_prime(q)) || continue
                for P ∈ E
                    (is_infinity(P) || P.y === nothing || iszero(P.y)) && continue
                    Q = scalar_mult(h, P, proposed)
                    is_infinity(Q) && continue
                    order = point_order(Q, proposed)
                    if q > best_order
                        best_order = q
                        best_curve = Montgomery{F}(b, a, q, h, (Int(Q.x), Int(Q.y)))
                        @info("New best $best_curve")
                    end
                end
            end
        end
    end
    isnothing(best_curve) ? @info("No suitable Montgomery curve with generator found.") :
                            @info("Returning best effort curve with subgroup order $best_order")
    return best_curve
end
