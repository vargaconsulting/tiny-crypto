using ..Macroes
using ..Utils

struct Weierstrass{T} <: Curve
    @define(T, π,a,b, order) # prime, a,b, order
    G::Point{T}  # Generator point (x, y)
end

function Base.show(io::IO, curve::Weierstrass{T}) where {T}
    @attach(curve, π, a, b, order, G)
    terms = []
    push!(terms, "x³")  # x³ is always present
    if a != 0
        push!(terms, "$a x")  # Include x term only if a ≠ 0
    end
    if b != 0
        push!(terms, "$b")  # Include constant term only if b ≠ 0
    end
    # Join the terms with " + " for a proper display
    curve_equation = join(terms, " + ")
    print(io, "Weierstrass curve: y² = $curve_equation mod $π with order: $order and generator point: $G")
end

# Point addition
function point_add(P::Point{T}, Q::Point{T}, curve::Weierstrass{T})::Point{T} where {T}
    if P.x === nothing
        return Q
    elseif Q.x === nothing
        return P
    elseif P.x == Q.x && P.y != Q.y
        return Point{T}(nothing, nothing)
    end

    π = curve.π
    λ = if P != Q
        mod((Q.y - P.y) * mod_inverse(Q.x - P.x, π), π)
    else
        num = mod(3 * powermod(P.x, 2, π) + curve.a, π)            ### FIXED
        den = mod_inverse(2 * P.y, π)
        mod(num * den, π)
    end

    x₃ = mod(powermod(λ, 2, π) - P.x - Q.x, π)                     ### FIXED
    y₃ = mod(λ * (P.x - x₃) - P.y, π)
    return Point{T}(x₃, y₃)
end


# Scalar multiplication
function scalar_mult(k::Integer, P::Point{T}, curve::Weierstrass{T}) where {T<:Integer}
    if k < 0
        return scalar_mult(-k, Point{T}(P.x, (-P.y) % curve.π), curve)
    end
    R, N = Point{T}(nothing, nothing), P
    while k > 0
        if k % 2 == 1
            R = point_add(R, N, curve)
        end
        N = point_add(N, N, curve)
        k ÷= 2
    end
    return R
end

function is_point_on_curve(P::Point{T}, curve::Weierstrass{T}) where {T}
    if P.x === nothing
        return true
    end
    π = curve.π
    lhs = mod(powermod(P.y, 2, π), π)                              ### FIXED
    rhs = mod(powermod(P.x, 3, π) + curve.a * P.x + curve.b, π)   ### FIXED
    return lhs == rhs
end


function curve_points(curve::Weierstrass{T}) where {T}
    @attach(curve, a,b,π,order)
    quadratic_residues = Set(powermod(y, 2, π) for y in 0:π-1)     ### FIXED
    points = []
    for x in 0:π-1
        y² = mod(powermod(x, 3, π) + a * x + b, π)                 ### FIXED
        if y² in quadratic_residues
            for y in 0:π-1
                if powermod(y, 2, π) == y²                         ### FIXED
                    push!(points, Point{T}(x, y))
                end
            end
        end
    end
    push!(points, Point{T}(nothing, nothing))
    return points
end


function is_generator(G::Point, curve::Weierstrass, n::Integer)
    T = typeof(curve.π)
    π = curve.π
    lhs = mod(powermod(G.y, 2, π), π)                              ### FIXED
    rhs = mod(powermod(G.x, 3, π) + curve.a * G.x + curve.b, π)   ### FIXED
    if lhs != rhs
        return false
    end
    return scalar_mult(convert(T, n), G, curve) == Point{T}(nothing, nothing)
end


function Weierstrass(Π::UnitRange{T}, A::UnitRange{T}, B::UnitRange{T} ) where {T}
    for π ∈ primes(Π)  # Generate primes within the range
        for a ∈ A  # Loop over curve parameter a
            for b ∈ B  # Loop over curve parameter b
                if mod(4a^3 + 27b^2, π) == 0 continue end  # Check for singularity
                proposed_curve = Weierstrass{T}(π, a, b, 0, Point{T}(nothing, nothing))
                E = curve_points(proposed_curve)  # Generate points on the curve
                N = length(E)  # Compute group order
                if !is_prime(N) continue end  # Skip if group order is not prime
                for P ∈ E
                    if P.x === nothing continue end  # Skip point at infinity
                    if is_generator(P, proposed_curve, N)
                        return Weierstrass{T}(π, a, b, N, Point{T}(P))
                    end
                end
            end
        end
    end
    error("No suitable Weierstrass curve found in given ranges.")
end

function Weierstrass{T}(πs::UnitRange, as::UnitRange, bs::UnitRange) where {T<:Integer}
    cast(x) = T(x)
    πs_cast = T(first(πs)):T(last(πs))
    as_cast = T(first(as)):T(last(as))
    bs_cast = T(first(bs)):T(last(bs))
    return Weierstrass(πs_cast, as_cast, bs_cast)
end
