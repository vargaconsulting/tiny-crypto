using TinyCrypto

const 𝔽₃₁ = Fp{UInt8, 31}

## Weierstrass curve
𝔽ₚ
curve = Weierstrass{𝔽₃₁}(6, 9, 37, 1, (0,3)) # direct definition of a curve from parameters on a prime field
# or directly: 
curve = Weierstrass(31, 6, 9, 37, 1, (0,3))  # ditto, different syntax
# Weierstrass curve: y² = x³ + 6x + 9 |𝔽₃₁ with order: 37 and 𝔾(0,3)


#or finding one    (𝔽₁-𝔽₂,   a₁-a₂, b₁-b₂)
curve = Weierstrass(30:300, 5:10, 4:10)
E = curve_points(curve)
S = subgroup_points(curve)



# prints out : 101-element Vector{Point{Fp{UInt128, 83}}}
# can it be changed to 101-element Vector{AffinePoint₁₂₈}} where 128 denotes bit width?
curve.G
O = infinity(curve)
O = identity(curve)



curve = Weierstrass(15:50, 5:10, 4:10)
E = curve_points(curve)
O = infinity(curve)
is_point_on_curve(E[1], curve)
is_infinity(O)


{(11,13),(14,4), ... ,(1,8),(16,16),(5,15),(1,9),(9,4),(13,5),(3,10)} ⊂ Weierstrass{𝔽₁₇}

∞
point_order(E[3], curve)
is_generator(E[1])

curve.point_order * curve.G

# Montgomery curve: 1y² = x³ + 5x² + x |𝔽₂₆₃ with order: 73 and 𝔾(4,197)

using TinyCrypto
curve = Montgomery(30:40, 8:40, 3:40) 

E = curve_points(curve)
O = infinity(curve)
curve.G
curve.order * curve.G
S = subgroup_points(curve)
[curve.order * P for P ∈ S]  # all ∞
# Montgomery curve: 2y² = x³ + 1x² + x |𝔽₁₀₁₃ with order: 241 and 𝔾(786,772)



using TinyCrypto
import TinyCrypto:is_complete, point_order, is_infinity, is_identity, identity
curve1 = TwistedEdwards(50:100, 1:20, 1:20)
curve2 = Edwards(50:100, 1:20)


curve = Edwards(83, 6, 11, 8, (3, 58))
curve = Edwards(20:40, 2:10) 


this is what it outputs: Edwards{Fp{UInt64, 0x0000000000000053}}(TwistedEdwards{𝔽₈₃}: 1x² + y² = 1 + 6x²y² | 𝔾(3,58), q = 11, h = 8, #E = 88)
# my expectation: Edwards{𝔽₈₃}: 1x² + y² = 1 + 6x²y² | 𝔾(3,58), q = 11, h = 8, #E = 88

is_square(curve.a, curve.π) # true
is_square(curve.d, curve.π) # false



using TinyCrypto
import TinyCrypto:is_complete, point_order, is_infinity, is_identity, identity

curve = TwistedEdwards(20:40, 2:10, 2:10) 
is_square(curve.a, curve.π) # true
is_square(curve.d, curve.π) # false
is_complete(curve)          # true

E = [AffinePoint(p, curve) for p ∈ curve_points(curve)]

subgroup_points(curve)
i = 1; @info "index: $i point $(E[i]) n * cofactor * $(E[i]) where n ∈ 0:#E(𝔽₃₇) - 1" ; [n*curve.cofactor * E[i] for n ∈ 0:curve.order-1]

## testing generator point
G = E[i]
subgroup = [n * G for n in 0:curve.order-1]
length(unique(subgroup)) == curve.order || error("Points not unique ⇒ G not generator.") # true
all(is_point_on_curve(P, curve) for P in subgroup) # false

for (i, P) in enumerate(subgroup)
    if !is_point_on_curve(P, curve)
        @warn "Bad point at n=$i: $P"
    end
end
F = 𝔽ₚ{UInt128, 37}
P = AffinePoint(Point{F}(21, 2), curve)
@info "is_point_on_curve(P, curve) where P=$P and $curve   $(is_point_on_curve(P, curve))"
@info "P + P = $(P+P) "

A = point_add(P.point, P.point, curve) 

@info "G = $i * $(E[i]) = $(i*E[i])"
p = E[4] 

# IDENTITY Check
I = AffinePoint(Point{F}(0, 1), curve)
is_infinity(I) # false 
is_identity(I) # true
∞ = AffinePoint(Point{F}(nothing, nothing), curve)
is_infinity(∞) # true
is_identity(∞) # false




[is_point_on_curve(P, curve) for P in subgroup]

## TODO: test if n* ∈ E
[ n*G  for n ∈ 0:curve.order-1]

point_order(E[i].point, curve)


O = Infinit2y(curve)
S = subgroup_points(curve)
G = curve.G
is_generator(G, curve, curve.order)

curve.order * curve.G

[curve.order * P for P ∈ S]  # all ∞
# Montgomery curve: 2y² = x³ + 1x² + x |𝔽₁₀₁₃ with order: 241 and 𝔾(786,772)


