using ..Macroes
using ..Utils
using ..Field: Fp, 𝔽ₚ, isprime

struct Weierstrass{F<:Fp} <: Curve
    @define(Int, π, a, b, order) # TODO: what If I just used integer model praramers?
    G::Point{F}
end

function Weierstrass{F}(a::Integer, b::Integer, order::Integer, G::Tuple{Integer,Integer}) where {F<:Fp}
    π =  Int(F.parameters[2])  # Just to extract the modulus
    curve = Weierstrass{F}(
        π, a, b, order, Point{F}(F(G[1]), F(G[2])) )
    return curve
end
function Weierstrass(π::T, a::Integer, b::Integer, order::Integer, G::Tuple{Integer,Integer}) where {T<:Integer}
    𝔽ₚ = Fp{unsigned(T), reinterpret(unsigned(T), π)}
    return Weierstrass{𝔽ₚ}(a,b,order, G)
end

function Base.show(io::IO, curve::Weierstrass{F}) where {F<:Fp}
    @attach(curve, π, a, b, order, G)
    @assert π == F.parameters[2]

    a_int, b_int = Int(a), Int(b)  # Extract plain integers

    # Format bold numbers using ANSI escape codes
    bold(x) = "\e[1m$(x)\e[22m"

    terms = ["x³"]
    a_int != 0 && push!(terms, "$(bold(a_int))x")
    b_int != 0 && push!(terms, bold(b_int))

    curve_eq = join(terms, " + ")

    P = F.parameters[2]
    subscript = join(Char(0x2080 + d) for d in reverse(digits(P)))
    field_str = "𝔽$subscript"

    gx = Int(G.x)
    gy = Int(G.y)
    G_str = "𝔾($(bold(gx)),$(bold(gy)))"

    print(io, "Weierstrass curve: y² = $curve_eq |$field_str with order: $(bold(order)) and $G_str")
end

# Weierstrass curve: y² = x³ + 1x + 30 |𝔽₃₁ with order: 31 and 𝔾(6,2)
# can I have the coefficients a 1,b 30 the order 31 printed bold?

function point_add(P::Point{F}, Q::Point{F}, curve::Weierstrass{F})::Point{F} where {F<:Fp}
    if is_infinity(P)
        return Q
    elseif is_infinity(Q)
        return P
    elseif P.x == Q.x && P.y != Q.y
        return infinity(F)
    end

    λ = if P != Q
        (Q.y - P.y) / (Q.x - P.x)
    else
        (F(3) * P.x^2 + curve.a) / (F(2) * P.y)
    end

    x3 = λ^2 - P.x - Q.x
    y3 = λ * (P.x - x3) - P.y
    return Point{F}(x3, y3)
end

function scalar_mult(k::Integer, P::Point{F}, curve::Weierstrass{F})::Point{F} where {F<:Fp}
    k < 0 && return scalar_mult(-k, Point{F}(P.x, -P.y), curve)
    R = infinity(F)
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

function is_point_on_curve(P::Point{F}, curve::Weierstrass{F}) where {F<:Fp}
    is_infinity(P) && return true
    lhs = P.y^2
    rhs = P.x^3 + curve.a * P.x + curve.b
    return lhs == rhs
end

function curve_points(curve::Weierstrass{F}) where {F<:Fp}
    @attach(curve, a, b, π, order) # << we have π here defined, why not just use it? 
    πᵢ = Int(F.parameters[2])

    points = Point{F}[]
    for x_raw in 0:πᵢ-1
        x = F(x_raw)
        rhs = x^3 + a * x + b

        found = false
        for y_raw in 0:πᵢ-1
            y = F(y_raw)
            if y^2 == rhs
                push!(points, Point{F}(x, y))
                found = true
            end
        end
    end

    push!(points, infinity(F))
    return points
end

function is_generator(G::Point{F}, curve::Weierstrass{F}, n::Integer) where {F<:Fp}
    is_point_on_curve(G, curve) || return false
    return scalar_mult(n, G, curve) == infinity(F)
end

function Weierstrass(Π::UnitRange, A::UnitRange, B::UnitRange, ::Val{T}=Val(UInt128)) where {T<:Unsigned}
    for π ∈ primes(Π)
        isprime(π) || continue
        for a ∈ A, b ∈ B
            is_singular = mod(4a^3 + 27b^2, π) == 0
            is_singular && continue
            F = Fp{T, π}
            proposed = Weierstrass{F}(a, b, 0, (0, 0))  # dummy order, dummy generator
            pts = curve_points(proposed)
            N = length(pts)
            
            isprime(N) || continue
            
            for P in pts
                if !is_infinity(P) && is_generator(P, proposed, N)
                    return Weierstrass{F}(a, b, N, (Int(P.x), Int(P.y)))
                end
            end
        end
    end
    @info("No suitable curve found.")
end

