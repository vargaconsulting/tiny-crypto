using Test
using TinyCrypto

@testset "Weierstrass curve sanity" begin
    curve = Weierstrass(97:103, 0:2, 6:100)
    G = curve.G

    @test is_point_on_curve(G, curve)
    @test is_generator(G, curve, curve.order)

    # 0 * G = ∞
    P = scalar_mult(0, G.point, curve)
    @test is_infinity(P) 
    P = 0 * G
    @test is_infinity(P) 

    # 1 * G = G
    P = scalar_mult(1, G.point, curve)
    @test P == G.point
    P = 1 * G
    @test P == G

    # 2 * G = G + G
    P₁ = scalar_mult(2, G.point, curve)
    P₂ = point_add(G.point, G.point, curve)
    @test P₁ == P₂
    P₁,P₂ = 2 * G, G + G
    @test P₁ == P₂

    # Scalar multiplication distributes over addition
    P₁ = scalar_mult(3, G.point, curve)
    P₂ = point_add(G.point, point_add(G.point, G.point, curve), curve)
    @test P₁ == P₂
    P₁, P₂ = 3 * G, G + G + G
    @test P₁ == P₂

    # Generator of prime order should cycle back to identity
    P = scalar_mult(curve.order, G.point, curve)
    @test is_infinity(P)
    P = curve.order * G
    @test is_infinity(P)
end
