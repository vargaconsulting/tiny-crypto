
using ..Macroes
using ..Utils
using ..Field: Fp, 𝔽ₚ, isprime

mutable struct Montgomery{F<:Fp} <: Curve
    @define(Int, π, A, B, order)
    G::ECPoint{F, Montgomery{F}}

    function Montgomery{F}(A::Integer, B::Integer, order::Integer, Gxy::Tuple{Integer,Integer}) where {F<:Fp}
        π = Int(F.parameters[2])
        self = new{F}(π, A, B, order)  # Incomplete init
        self.G = ECPoint(Point{F}(F(Gxy[1]), F(Gxy[2])), self)  # Complete init
        return self
    end
end
function Montgomery(π::T, A::Integer, B::Integer, order::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    𝔽ₚ = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return Montgomery{𝔽ₚ}(A, B, order, G)
end
curve_equation_rhs(P::Point{F}, curve::Montgomery{F}) where {F<:Fp} = P.x^3 + curve.A * P.x^2 + P.x
curve_equation_lhs(P::Point{F}, curve::Montgomery{F}) where {F<:Fp} = curve.B * P.y^2
function is_singular(curve::Montgomery{F}) where {F<:Fp}
    @attach(curve, A, B, π)
    return B == 0 || mod(A^2 - 4, π) == 0
end

function Base.show(io::IO, curve::Montgomery{F}) where {F<:Fp}
    @attach(curve, π, A, B, order, G)
    A_int, B_int = Int(A), Int(B)
    bold(x) = "\e[1m$(x)\e[22m"

    terms = ["x³"]
    A_int != 0 && push!(terms, "$(bold(A_int))x²")
    push!(terms, "x")
    rhs_str = join(terms, " + ")
    lhs_str = "$(bold(B_int))y²"

    curve_eq = "$lhs_str = $rhs_str"

    subscript = join(Char(0x2080 + d) for d in reverse(digits(π)))
    field_str = "𝔽$subscript"

    gx = Int(G.point.x)
    gy = Int(G.point.y)
    G_str = "𝔾($(bold(gx)),$(bold(gy)))"

    print(io, "Montgomery curve: $curve_eq |$field_str with order: $(bold(order)) and $G_str")
end

function point_add(P::Point{F}, Q::Point{F}, curve::Montgomery{F})::Point{F} where {F<:Fp}
    if is_infinity(P)
        return Q
    elseif is_infinity(Q)
        return P
    elseif P.x == Q.x
        if P.y + Q.y == F(0)  # additive inverse
            return infinity(F)
        elseif P == Q && P.y == F(0)  # vertical tangent in doubling
            return infinity(F)
        end
    end

    λ = if P != Q
        (Q.y - P.y) / (Q.x - P.x)
    else
        # Doubling
        num = F(3) * P.x^2 + F(2) * curve.A * P.x + F(1)
        den = F(2) * curve.B * P.y
        num / den
    end

    x3 = curve.B * λ^2 - curve.A - P.x - Q.x
    y3 = λ * (P.x - x3) - P.y
    return Point{F}(x3, y3)
end

function Montgomery(Π::UnitRange, A::UnitRange, B::UnitRange, max_cofactor::Int = 8, ::Val{T}=Val(UInt128)) where {T<:Unsigned}
    best_curve, best_order = nothing, 0
    for π ∈ primes(Π)
        for a in filter(<(π), A), b in filter(<(π), B)
            F = 𝔽ₚ{T, π}
            proposed = Montgomery{F}(a, b, 0, (0, 0))
            is_singular(proposed) && continue
            E = curve_points(proposed)
            N = length(E)
            for P ∈ E
                (P.y === nothing || iszero(P.y)) && continue
                for q in reverse(sort(prime_factors(N)))
                    h = div(N, q)
                    h > max_cofactor && continue
                    G = h * ECPoint(P, proposed)
                    (is_infinity(G) || iszero(G.point.y)) && continue
                    is_point_on_curve(G, proposed) || continue
                    curve = Montgomery{F}(a, b, q, (Int(G.point.x), Int(G.point.y)))
                    if q > best_order
                        best_order = q
                        best_curve = curve
                        @info("New best: subgroup order $q (cofactor $h) $curve")
                    end
                end
            end
        end
    end
    isnothing(best_curve) ? @info("No suitable Montgomery curve with generator found.") :
                            @info("Returning best effort curve with subgroup order $best_order")
    return best_curve
end
