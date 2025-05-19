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


function fiat_shamir(t::Point{Int}, y::Point{Int}, q::Int)
    bytes = join([t.x, t.y, y.x, y.y])
    hash = sha8(bytes)
    return Int(hash % q)
end
function fiat_shamir(t::Point{Int}, y::Point{Int}, q::Int, msg::AbstractString)
    bytes = join([t.x, t.y, y.x, y.y, msg])
    h = sha8(bytes)
    return Int(h % q)
end

function schnorr_sign(curve::Weierstrass{Int}, x::Int)
    @attach(curve, order, G)
    y = scalar_mult(x, G, curve)

    nonce = rand(1:order-1)
    t = scalar_mult(nonce, G, curve)
    c = fiat_shamir(t, y, order)
    s = mod(nonce + c * x, order)

    return t, s, y
end

function schnorr_verify(curve::Weierstrass{Int}, y::Point{Int}, t::Point{Int}, s::Int)
    @attach(curve, order, G)
    c = fiat_shamir(t, y, order)
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

# Non-Interactive Schnorr ZKPoK (via Fiat–Shamir)
x = rand(1:curve.order-1)                 # Prover's secret key: random scalar in the group
t, s, y = schnorr_sign(curve, x)          # Prover computes commitment t = G^r, Fiat–Shamir challenge c = H(t,y), response s = r + c·x
valid = schnorr_verify(curve, y, t, s)    # Verifier recomputes c = H(t,y), checks that G^s == t * y^c
