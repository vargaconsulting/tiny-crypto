using Test
using TinyCrypto
import TinyCrypto: is_identity, is_infinity, is_singular, inverse,
                   point_neg, point_add, curve_points, scalar_mult,
                   is_point_on_curve, is_generator, subgroup_points

@testset "Weierstrass Curve Test Suite" begin
    # Example Weierstrass curve: y² = x³ + 2x + 2 over 𝔽₁₇ with generator (5,1)
    π, a, b, q, h = 17, 2, 2, 19, 1
    Gxy = (5, 1)
    curve = Weierstrass(π, a, b, q, h, Gxy)
    F = typeof(curve.G.point.x)

    Gᶜ = curve.G
    G = Gᶜ.point
    ∞ = infinity(curve)
    ∞ᶜ = ECPoint(∞, curve)
    S = subgroup_points(curve)

    @test is_infinity(5 * ∞ᶜ)

    @testset "Sanity checks and scalar multiplication" begin
        @test is_point_on_curve(Gᶜ, curve)
        @test is_generator(G, curve, curve.order)

        @test is_identity(scalar_mult(0, G, curve), curve)
        @test is_identity(0 * Gᶜ)

        @test scalar_mult(1, G, curve) == G
        @test 1 * Gᶜ == Gᶜ

        @test scalar_mult(2, G, curve) == point_add(G, G, curve)
        @test 2 * Gᶜ == Gᶜ + Gᶜ

        @test is_identity(curve.order * Gᶜ)
    end

    @testset "Subgroup matches [n]Gᶜ for n ∈ 0:q-1" begin
        expected = Set(n * Gᶜ for n in 0:curve.order - 1)
        actual   = Set(S)
        @test expected == actual
    end

    @testset "Point + identity = point" begin
        @test point_add(G, ∞, curve) == G
        @test point_add(∞, G, curve) == G
    end

    @testset "Point + inverse = identity" begin
        N = point_neg(G, curve)
        @test is_identity(point_add(G, N, curve), curve)
    end

    @testset "Inverse properties" begin
        for Pᶜ in S
            Qᶜ = inverse(Pᶜ)
            @test is_identity(Pᶜ + Qᶜ)
            @test inverse(inverse(Pᶜ)) == Pᶜ
            if !is_infinity(Pᶜ)
                @test Qᶜ.point.x == Pᶜ.point.x
                @test Qᶜ.point.y == -Pᶜ.point.y
            end
        end
    end

    @testset "Group law: [m]G + [n]G == [(m+n)%q]G" begin
        for m in 0:q-1, n in 0:q-1
            @test m * Gᶜ + n * Gᶜ == mod(m+n, q) * Gᶜ
        end
    end

    @testset "All curve points lie on the curve" begin
        for P in curve_points(curve)
            @test is_point_on_curve(ECPoint(P, curve), curve)
        end
    end

    @testset "Subgroup checks" begin
        @test all(is_identity, [q * Pᶜ for Pᶜ in S])
        @test length(S) == q
        @test length(unique(S)) == q
    end

    @testset "Reject invalid points" begin
        bad = Point{F}(F(42), F(1337))
        @test !is_point_on_curve(ECPoint(bad, curve), curve)
    end

    @testset "Doubling consistency" begin
        D₁ = point_add(G, G, curve)
        D₂ = scalar_mult(2, G, curve)
        @test D₁ == D₂
    end

    @testset "Scalar multiplication with ∞" begin
        for k in 0:10
            @test is_infinity(k * ∞ᶜ)
        end
    end

    @testset "Curve is not singular" begin
        @test !is_singular(curve)
    end

    @testset "Negation operator -Pᶜ" begin
        for Pᶜ in S
            @test -Pᶜ == inverse(Pᶜ)
        end
    end

    @testset "Negation specific cases" begin
        Q = point_neg(∞, curve)
        @test Q == ∞

        G⁻ = point_neg(G, curve)
        @test is_identity(point_add(G, G⁻, curve), curve)

        G⁻⁻ = point_neg(G⁻, curve)
        @test G⁻⁻ == G
    end

    @testset "Associativity sanity check" begin
        Dᶜ = ECPoint(point_add(G, G, curve), curve)
        G⁻ = point_neg(G, curve)
        A = point_add(Dᶜ.point, G⁻, curve)
        B = point_add(G, point_add(G, G⁻, curve), curve)
        @test A == B
    end

    @testset "Infinity propagation" begin
        R₁ = point_add(G, ∞, curve)
        R₂ = point_add(∞, G, curve)
        @test R₁ == G
        @test R₂ == G
    end
end
