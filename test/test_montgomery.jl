using Test
using TinyCrypto
import TinyCrypto: is_identity, is_infinity, inverse, is_singular, point_neg, point_add

@testset "Montgomery Curve Test Suite" begin
    # Shared curve setup Montgomery{𝔽₃₇}: 8y² = x³ + 3x² + x | 𝔾(15,2), q = 11, h = 4, #E = 44
    π, B, A, order, cofactor = 37, 8, 3, 11, 4
    Gˣʸ = (15, 2)
    curve = Montgomery(π, B, A, order, cofactor, Gˣʸ)

    Gᶜ = curve.G
    G   = Gᶜ.point
    F   = typeof(G.x)

    I  = identity(curve).point
    Iᶜ = AffinePoint(infinity(curve), curve)
    ∞ᶜ = AffinePoint(infinity(curve), curve)
    Gᶜ + Gᶜ + Gᶜ + Gᶜ
    4Gᶜ
    curve
    curve_points(curve)
    S   = subgroup_points(curve)

    @testset "Subgroup matches [n]Gᶜ for n ∈ 0:q-1" begin
        expected = Set(n * Gᶜ for n in 0:curve.order - 1)
        actual   = Set(S)
        @test expected == actual
    end

    @testset "Sanity checks and scalar multiplication" begin
        @test is_point_on_curve(Gᶜ, curve)
        @test is_generator(G, curve, curve.order)

        @test is_infinity(scalar_mult(0, G, curve))
        @test is_infinity(0 * Gᶜ)

        @test scalar_mult(1, G, curve) == G
        @test 1 * Gᶜ == Gᶜ

        @test scalar_mult(2, G, curve) == point_add(G, G, curve)
        @test 2 * Gᶜ == Gᶜ + Gᶜ

        @test scalar_mult(3, G, curve) == (3 * Gᶜ).point
        @test is_infinity(curve.order * Gᶜ)

        if curve.order > 5
            @test is_point_on_curve(2 * Gᶜ, curve)
            @test is_point_on_curve(3 * Gᶜ, curve)
            @test 3 * Gᶜ == 2 * Gᶜ + Gᶜ
        end
    end

    @testset "Point + identity = point" begin # OK
        R₁ = point_add(G, I, curve)
        R₂ = point_add(I, G, curve)
        @test R₁ == G
        @test R₂ == G
    end

    @testset "Point + inverse = identity" begin
        Pᶜ = Gᶜ
        Nᶜ = AffinePoint(point_neg(Pᶜ.point, curve), curve)
        Rᶜ = AffinePoint(point_add(Pᶜ.point, Nᶜ.point, curve), curve)
        @test is_infinity(Rᶜ)
    end

    @testset "Inverse properties" begin
        @test is_infinity(inverse(∞ᶜ))
        for Pᶜ in S
            Qᶜ = inverse(Pᶜ)
            @test is_infinity(Pᶜ + Qᶜ)
            @test inverse(inverse(Pᶜ)) == Pᶜ
            if !is_infinity(Pᶜ)
                @test Qᶜ.point.x == Pᶜ.point.x
                @test Qᶜ.point.y == -Pᶜ.point.y
            end
        end
    end

    @testset "All curve points satisfy the curve equation" begin
        for P in curve_points(curve)
            @test is_point_on_curve(AffinePoint(P, curve), curve)
        end
    end

    @testset "Scalar multiplication: [n]G ∈ subgroup, ∀ n < order" begin
        for n in 0:curve.order - 1
            @test is_point_on_curve(n * Gᶜ, curve)
        end
    end

    @testset "Montgomery group law: [m]G + [n]G == [(m+n)%r]G" begin
        r = curve.order
        for m in 0:r-1, n in 0:r-1
            P = m * Gᶜ
            Q = n * Gᶜ
            @test P + Q == mod(m + n, r) * Gᶜ
        end
    end

    @testset "Curve singularity check" begin
        @test !is_singular(curve)
    end

    @testset "Reject invalid points" begin
        bad = Point{F}(F(123), F(456))
        @test !is_point_on_curve(AffinePoint(bad, curve), curve)
    end

    @testset "Doubling edge case: vertical tangent" begin
        for P in curve_points(curve)
            if P.y !== nothing && P.y == F(0)
                @test is_infinity(point_add(P, P, curve))
            end
        end
    end

    @testset "Doubling matches scalar multiplication" begin
        @test point_add(G, G, curve) == scalar_mult(2, G, curve)
    end

    @testset "Scalar mult with ∞" begin
        for n in 0:10
            @test is_infinity(n * ∞ᶜ)
        end
    end

    @testset "Operator overload for negation" begin
        for P in S
            @test -P == inverse(P)
        end
    end
end
