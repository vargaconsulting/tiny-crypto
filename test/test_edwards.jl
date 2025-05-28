using Test
using TinyCrypto
import TinyCrypto: is_identity, is_infinity, inverse, is_singular, point_neg, point_add

@testset "Edwards Curve Test Suite" begin
    # Shared curve setup Twisted Edwards curve: 1x² + y² = 1 + 15x²y² |𝔽₃₇ with order: 11 cofactor: 4 and 𝔾(18,17)
    π, a, d, order, cofactor = 37, 1, 15, 11, 4
    Gˣʸ    = (18, 17)   
    curve = TwistedEdwards(π, a, d, order, cofactor, Gˣʸ)
    F = typeof(curve.G.point.x)

    Gᶜ = curve.G
    G = Gᶜ.point
    I = Point{F}(F(0), F(1))
    Iᶜ = ECPoint(I, curve)
    ∞ = Point{F}(nothing, nothing)
    ∞ᶜ = ECPoint(∞, curve)
    S = subgroup_points(curve)

    @test is_infinity(5 * ∞ᶜ)

    @testset "Subgroup matches [n]Gᶜ for n ∈ 0:q-1" begin
        expected = Set(n * Gᶜ for n in 0:curve.order - 1)
        actual   = Set(S)
        @test expected == actual
    end

    @testset "Sanity checks and scalar multiplication" begin
        @test is_point_on_curve(Gᶜ, curve)
        @test is_generator(G, curve, curve.order)

        @test is_identity(scalar_mult(0, G, curve), curve)
        @test is_identity(0 * Gᶜ)

        @test scalar_mult(1, G, curve) == G
        @test 1 * Gᶜ == Gᶜ

        @test scalar_mult(2, G, curve) == point_add(G, G, curve)
        @test 2 * Gᶜ == Gᶜ + Gᶜ

        @test scalar_mult(3, G, curve) == (3 * Gᶜ).point

        @test is_identity(curve.order * Gᶜ)
    end

    @testset "Point + identity = point" begin
        Pᶜ = Gᶜ
        R₁ = point_add(Pᶜ.point, Iᶜ.point, curve)
        R₂ = point_add(Iᶜ.point, Pᶜ.point, curve)
        @test R₁ == G
        @test R₂ == G
    end

    @testset "Point + inverse = identity" begin
        Pᶜ = Gᶜ
        Nᶜ = ECPoint(point_neg(Pᶜ.point, curve), curve)
        R₁ᶜ = ECPoint(point_add(Pᶜ.point, Nᶜ.point, curve), curve)
        @test is_identity(R₁ᶜ)
    end

    @testset "Inverse properties" begin
        @test is_infinity(inverse(∞ᶜ))
        for Pᶜ in S
            Qᶜ = inverse(Pᶜ)
            @test is_identity(Pᶜ + Qᶜ)
            @test inverse(inverse(Pᶜ)) == Pᶜ
            if !is_infinity(Pᶜ)
                @test Qᶜ.point.x == -Pᶜ.point.x
                @test Qᶜ.point.y == Pᶜ.point.y
            end
        end
    end

    @testset "Group law: [m]G + [n]G == [(m+n)%r]G" begin
        r = curve.order
        for m in 0:r-1, n in 0:r-1
            Pᶜ = m * Gᶜ
            Qᶜ = n * Gᶜ
            @test Pᶜ + Qᶜ == mod(m+n, r) * Gᶜ
        end
    end

    @testset "All curve points lie on the curve" begin
        for P in curve_points(curve)
            @test is_point_on_curve(ECPoint(P, curve), curve)
        end
    end

    @testset "Subgroup checks" begin
        @test all(is_identity, [curve.order * Pᶜ for Pᶜ in S])
        @test length(S) == curve.order
        @test length(unique(S)) == curve.order
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
        # Identity negation
        Q = point_neg(I, curve)
        @test Q == I

        # Infinity negation
        Q = point_neg(∞, curve)
        @test Q == ∞

        # Generator negation
        G⁻ = point_neg(G, curve)
        @test is_identity(ECPoint(point_add(G, G⁻, curve), curve))

        # Double negation
        G⁻⁻ = point_neg(G⁻, curve)
        @test G⁻⁻ == G

        # Arbitrary valid point
        P = Point{F}(F(5), F(7))
        if is_point_on_curve(P, curve)
            Pᶜ = ECPoint(P, curve)
            P⁻ᶜ = ECPoint(point_neg(P, curve), curve)
            sumᶜ = ECPoint(point_add(P, P⁻ᶜ.point, curve), curve)
            @test is_identity(sumᶜ)
        end
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

