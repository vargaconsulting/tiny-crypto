
using TinyCrypto

# ─────────────────────────────────────────────────────────────────────────────
# ECDSA Signature and Public Key Recovery Example
# ─────────────────────────────────────────────────────────────────────────────

priv, pub = genkey(curve)                     # Key generation
msg = "hello ethereum"                        # Message to sign
signature = sign(curve, priv, msg)            # returns a NamedTuple (r = ..., s = ..., v = ...)
is_valid = verify(curve, pub, signature, msg) # Signature verification true
recovered = ecrecover(curve, msg, signature)  # Public key recovery from signature
@assert pub == recovered                      # Check if recovered public key matches originald

# ─────────────────────────────────────────────────────────────────────────────
# The following steps (0–7) manually walk through the ECDSA signing and verification 
# process to illustrate what happens under the hood. The goal is to help understand 
# the mathematical structure behind the cryptographic operations shown above.
# ─────────────────────────────────────────────────────────────────────────────

# 1. Define a toy Weierstrass curve over 𝔽₃₁ with G of order 37
curve = Weierstrass(31, 6, 9, 37, 1, (0, 3))
# Weierstrass{𝔽₃₁}: y² = x³ + 6x + 9 | G = (0,3), order = 37, cofactor = 1

# 2. Create a mock blockchain transaction
tx = "hello ethereum"

# 3. Hash it (you fake Keccak256 with H₈ → SHA256 first 8 bytes)
z = Int(H₈(tx)) % curve.order  # H₈ returns UInt8; convert to Int for modular arithmetic

# 4. Generate keypair
priv = rand(1:curve.order - 1)          # private key is a uniformly random integer in [1, curve order − 1]
pub  = priv * curve.G                   # public key is the scalar multiple of the generator: pub = priv · G

r,s,v = 0,0,0
# 5. Sign: create ephemeral scalar `k` and calculate signature (r, s)
while true
    k = rand(1:curve.order - 1)             # draw a random ephemeral scalar used once per signature
    R = k * curve.G                         # then construct an ephemeral public point R = k·G on the curve
    v = isodd(R.point.y.val) ? 1 : 0        # and compute the parity (Ethereum) bit from LSB of y which is used for ECrecover

    r = mod(R.point.x, curve.order)         # use x-coordinate of R as first signature component
    k_inv = invmod(k, curve.order)          # compute modular inverse of ephemeral scalar
    s = mod(k_inv * (z + r * priv), curve.order)  # second signature component s = k⁻¹(z + r·priv) mod n

    if r ≠ 0 && s ≠ 0 break end             # loop until suitable `k` is found
end

# 6. The signature is (r, s, v + 27)   `27` roots back to early Bitcoin compact signatures, parity bit is used on Ethereum
# whereas Bitcoin uses recovery ID : = {0,1,2,3}  
println("Signature: (r = $r, s = $s, v = $(v + 27))")

# 7. Verification
s_inv = invmod(s, curve.order)               # compute inverse of signature scalar s⁻¹ mod n
u₁ = z * s_inv % curve.order                 # u₁ = z · s⁻¹ mod n
u₂ = r * s_inv % curve.order                 # u₂ = r · s⁻¹ mod n
P = u₁ * curve.G + u₂ * pub                  # recover point P = u₁·G + u₂·Q (should equal R)

is_valid = mod(P.point.x, curve.order) == r  # check if recovered R.x matches r component
println("Signature valid? ", is_valid)
