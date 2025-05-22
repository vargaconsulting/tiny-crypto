using TinyCrypto
using Test

const F31 = Fp{UInt8, 31}

@testset "Fp finite field" begin
    a = F31(7)
    b = F31(13)

    # Basic operations
    @test a + b == F31(20)
    @test a - b == F31(25)
    @test a * b == F31(29)
    @test F31(5) * F31(17) == F31(23)
    @test F31(3) * F31(4) == F31(12)
    @test inv(a) * a == one(F31)
    @test a ^ 0 == one(F31)
    @test a ^ 3 == a * a * a

    # Constructor checks
    @test_throws ArgumentError Fp{UInt8, 30}(1)  # Not a prime
    @test F31(31) == F31(0)                      # wrap
    @test F31(32) == F31(1)                      # wrap

    # Inversion edge case
    @test_throws DivideError inv(F31(0))

    # Shift operators (bitwise, modulo-safe)
    @test F31(1) << 1 == F31(2)
    @test F31(2) << 2 == F31(8)
    @test F31(4) >> 1 == F31(2)
    @test F31(3) << 5 == F31(3 << 5 % 31)  # large shift wraparound

    # Identities and properties
    @test zero(F31) == F31(0)
    @test one(F31) == F31(1)
    @test iszero(F31(0))
    @test isone(F31(1))
    @test -F31(5) == F31(26)
    @test abs(F31(18)) == F31(18)

    # Powers
    @test F31(2)^5 == F31(32) == F31(1)  # 32 % 31 == 1
    @test F31(3)^30 == one(F31)         # Fermat's little theorem: a^(p-1) ≡ 1 mod p
end
