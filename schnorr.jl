include("curve.jl")
include("sha.jl")
using .SHA8: sha8

function schnorr_commit(curve::Weierstrass{Int}, x::Int)
    @attach(curve, order, G)
    nonce = rand(1:order-1)
    t = scalar_mult(nonce, G, curve)
    return t, nonce
end

function schnorr_response(curve::Weierstrass{Int}, x::Int, nonce::Int, c::Int)
    @attach(curve, order)
    return mod(nonce + c * x, order)
end

function schnorr_verify(curve::Weierstrass{Int}, y::Point{Int}, t::Point{Int}, c::Int, s::Int)
    @attach(curve, order, G)
    lhs = scalar_mult(s, G, curve)
    rhs = point_add(t, scalar_mult(c, y, curve), curve)
    return lhs == rhs
end

# Interactive Schnorr Zero-Knowledge Proof of Knowledge
curve = Weierstrass(80:100, 1:5, 1:5)       # Generate a random prime field Weierstrass curve with small parameters
x = rand(1:curve.order-1)                   # Prover's secret: a random scalar in the group
y = scalar_mult(x, curve.G, curve)          # Public key y = G^x, i.e., scalar multiplication of generator G by x

t, nonce = schnorr_commit(curve, x)         # Prover generates commitment t = G^r, using a random nonce r
c = rand(1:curve.order-1)                   # Verifier generates random challenge c (in practice, could be a Fiat–Shamir hash)
s = schnorr_response(curve, x, nonce, c)    # Prover computes response s = r + cx mod order (modular arithmetic)

is_valid = schnorr_verify(curve, y, t, c, s)  # Verifier checks that G^s == t * y^c, confirming knowledge of x without revealing it
