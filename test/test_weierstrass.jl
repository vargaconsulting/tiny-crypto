using Test
using TinyCrypto

@testset "Weierstrass curve sanity" begin
    curve = Weierstrass(97:103, 0:2, 6:100)
    G = curve.G

    @test is_point_on_curve(G, curve)
    @test is_generator(G, curve, curve.order)

    # 0 * G = ∞
    P0 = scalar_mult(0, G, curve)
    @test P0 == Point{Int}(nothing, nothing)

    # 1 * G = G
    P1 = scalar_mult(1, G, curve)
    @test P1 == G

    # 2 * G = G + G
    P2a = scalar_mult(2, G, curve)
    P2b = point_add(G, G, curve)
    @test P2a == P2b

    # Scalar multiplication distributes over addition
    P3a = scalar_mult(3, G, curve)
    P3b = point_add(G, point_add(G, G, curve), curve)
    @test P3a == P3b

    # Generator of prime order should cycle back to identity
    Pn = scalar_mult(curve.order, G, curve)
    @test Pn == Point{Int}(nothing, nothing)
end
