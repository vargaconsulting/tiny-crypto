#   Copyright © <2010-2025> Varga Consulting, Toronto, On     info@vargaconsulting.ca
macro attach(st, fields...)
	block = Expr(:block)
	for f in fields
		push!(block.args, :($f = $st.$f) )
    end
	esc(:($block))
end
macro define(T, fields...)
	block = Expr(:block)
	for f in fields
		push!(block.args, :($f ::$T ))
    end
	return esc(:($block))
end

function is_prime(n::Int)::Bool
    if n < 2 return false end
    for i in 2:floor(Int, sqrt(n))
        if n % i == 0
            return false
        end
    end
    return true
end

function primes(range::UnitRange{Int})
    return [n for n in range if is_prime(n)]
end

function is_quadratic_residue(a::Int, p::Int)::Bool
    if a % p == 0
        return true  # 0 is trivially a quadratic residue
    end
    result = powermod(a, (p - 1) ÷ 2, p)
    return result == 1
end
function mod_inverse(a::T, p::T) where T # Modular inverse
    return powermod(a, p - 2, p)         # Fermat's little theorem
end

struct Point{T} # Parametric type for Point
    x::Union{T, Nothing}
    y::Union{T, Nothing}
end
function Base.show(io::IO, P::Point{T}) where T
    if P.x === nothing || P.y === nothing
        print(io, "(∞,∞)")
    else
        print(io, "(", P.x, ",", P.y, ")")
    end
end

abstract type Curve end
struct Weierstrass{T} <: Curve
    @define(T, π,a,b, order) # prime, a,b, order
    G::Point{T}  # Generator point (x, y)
end

function Base.show(io::IO, curve::Weierstrass{T}) where T
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
function point_add(P::Point{Int}, Q::Point{Int}, curve::Weierstrass{Int})::Point{Int}
    if P.x === nothing  # P is the point at infinity
        return Q
    elseif Q.x === nothing  # Q is the point at infinity
        return P
    elseif P.x == Q.x && P.y != Q.y  # Vertical line
        return Point{Int}(nothing, nothing)  # Point at infinity
    end

    λ = if P != Q  # Case 1: Adding distinct points
        mod((Q.y - P.y) * mod_inverse(Q.x - P.x, curve.π), curve.π)
    else  # Case 2: Doubling the point
        mod((3 * P.x^2 + curve.a) * mod_inverse(2 * P.y, curve.π), curve.π)
    end

    x₃ = mod(λ^2 - P.x - Q.x, curve.π)
    y₃ = mod(λ * (P.x - x₃) - P.y, curve.π)
    return Point{Int}(x₃, y₃)
end

# Scalar multiplication
function scalar_mult(k::Int, P::Point{T}, curve::Weierstrass{T}) where {T}
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
    if P.x === nothing  # Point at infinity is always valid
        return true
    end
    lhs = (P.y^2) % curve.π
    rhs = (P.x^3 + curve.a * P.x + curve.b) % curve.π
    return lhs == rhs
end

function curve_points(curve::Weierstrass{T}) where {T}
    @attach(curve, a,b,π,order)
    quadratic_residues, points = Set(mod(y^2, π) for y in 0:π-1), []
    for x in 0:π-1 # Compute y^2 = x^3 + ax + b mod p
        y² = (x^3 + a * x + b) % π      # compute y² and check if
        if y² in quadratic_residues
            for y in 0:π-1 # Check if y^2 is a quadratic residue (has a modular square root)
                if (y^2 % π) == y²
                    push!(points, Point(x, y))  # Add (x, y)
                end
            end
        end
    end
    push!(points, Point{T}(nothing,nothing)) # Add the point at infinity
    return points
end

function is_generator(G::Point{T}, curve::Weierstrass{T}, n::T) where {T}
    # Check if G satisfies the curve equation
    if (G.y^2 % curve.π) != (G.x^3 + curve.a * G.x + curve.b) % curve.π
        return false # Point is not on the curve
    end
    # Check if n * G = point at infinity
    if scalar_mult(n, G, curve) != Point{T}(nothing, nothing)
        return false # point does not have order n
    end
    return true # valid generator point
end

function Weierstrass(Π::UnitRange{T}, A::UnitRange{T}, B::UnitRange{T} ) where T
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
                        return Weierstrass{T}(π, a, b, N, P)
                    end
                end
            end
        end
    end
end

## Edwards curve: x² + y² = 1 + dx²y²
struct Edwards{T} <: Curve
    @define(T, π, d, order)
    G::Point{T}
end

function Base.show(io::IO, curve::Edwards{T}) where T
    @attach(curve, π, d, order, G)
    print(io, "Edwards curve: x^2 + y^2 = 1 + $d x^2 y^2 mod $π with order: $order and generator point: $G")
end

function point_add(P::Point{T}, Q::Point{T}, curve::Edwards{T}) where T
    @attach(curve, π, d)
    if P.x === nothing || P.y === nothing
        return Q
    elseif Q.x === nothing || Q.y === nothing
        return P
    end

    x1, y1 = P.x, P.y
    x2, y2 = Q.x, Q.y

    denom_x = mod(1 + d * x1 * x2 * y1 * y2, π)
    if denom_x == 0 return Point{T}(nothing, nothing) end
    x3 = mod((x1 * y2 + y1 * x2) * mod_inverse(denom_x, π), π)

    denom_y = mod(1 - d * x1 * x2 * y1 * y2, π)
    if denom_y == 0 return Point{T}(nothing, nothing) end
    y3 = mod((y1 * y2 - x1 * x2) * mod_inverse(denom_y, π), π)

    return Point{T}(x3, y3)
end


function scalar_mult(k::Int, P::Point{T}, curve::Edwards{T}) where T
    R, N = Point{T}(0, 1), P  # identity = (0, 1)
    while k > 0
        if k & 1 == 1
            R = point_add(R, N, curve)
        end
        N = point_add(N, N, curve)
        k >>>= 1
    end
    return R
end

function curve_points(curve::Edwards{T}) where {T}
    @attach(curve, π, d)
    points = Point{T}[]
    for x in 0:π-1
        x² = mod(x^2, π)
        for y in 0:π-1
            y² = mod(y^2, π)
            lhs = mod(x² + y², π)
            rhs = mod(1 + d * x² * y², π)
            if lhs == rhs
                push!(points, Point(x, y))
            end
        end
    end
    return points
end
function is_point_on_curve(P::Point{T}, curve::Edwards{T}) where {T}
    if P.x === nothing || P.y === nothing
        return false  # Edwards curves don't include the point at infinity
    end
    @attach(curve, π, d)
    x² = mod(P.x^2, π)
    y² = mod(P.y^2, π)
    lhs = mod(x² + y², π)
    rhs = mod(1 + d * x² * y², π)
    return lhs == rhs
end

function is_generator(G::Point{T}, curve::Edwards{T}, n::T) where {T}
    # Check if the point lies on the curve
    if !is_point_on_curve(G, curve)
        return false
    end
    # Check if [n]G == identity element (0, 1)
    return scalar_mult(n, G, curve) == Point{T}(0, 1)
end

function Edwards(Π::UnitRange{T}, D::UnitRange{T}) where T
    for π ∈ primes(Π)
        for d ∈ D
            if d == 0 || d == 1 continue end
            proposed_curve = Edwards{T}(π, d, 0, Point{T}(0, 1))  # Start with identity
            E = curve_points(proposed_curve)  # Generate points on the curve
            N = length(E)  # Compute group order
            for h in 1:8
                r, rem = divrem(N, h)
                if rem != 0 || !is_prime(r)
                    continue
                end
                for P in shuffle(E)
                    if P.x === nothing || P.y === nothing continue end
                    if is_generator(P, proposed_curve, r)
                        return Edwards{T}(π, d, r, P)
                    end
                end
                
            end
        end
    end
    println("No Edwards curve with prime group order found")
end

