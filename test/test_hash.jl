using Test
using TinyCrypto: H, H8, H16, H₈, H₁₆

@testset "Hash - TinyHash API" begin

    @testset "Determinism" begin
        @test H("msg") == H("msg")
        @test H8("abc") == H8("abc")
        @test H16("abc") == H16("abc")
    end

    @testset "Types & Output Range" begin
        @test H8("") isa UInt8
        @test H16("") isa UInt16
        @test 0x00 <= H8("abc") <= 0xff
        @test 0x0000 <= H16("abc") <= 0xffff
    end

    @testset "Unicode and Multibyte Strings" begin
        @test H("𝒞rypto") isa UInt8
        @test H16("éçüΩ") != H16("ecuo")
    end

    @testset "Empty input" begin
        @test H8("") == H(collect(UInt8, ""))
        @test H16("") == H16(UInt8[])
    end

    @testset "Single and Repeated Characters" begin
        @test H8("a") == H("a")
        @test H8("aaaaaaaa") != H8("aaaaaaaaa")  # should trigger padding into another block
        @test H16("z"^1000) isa UInt16
    end

    @testset "Vectors of UInt8" begin
        @test H8(UInt8[0x61, 0x62, 0x63]) == H8("abc")
        @test H16(UInt8[0x61, 0x62, 0x63]) == H16("abc")
    end

    @testset "Aliases" begin
        str = "hello world"
        @test H(str) == H8(str)
        @test H₈(str) == H8(str)
        @test H₁₆(str) == H16(str)
    end

    @testset "Corner padding boundaries" begin
        # Pad to block sizes: exact 7, 8, 14, 15, etc.
        for n in (0, 1, 6, 7, 8, 13, 14, 15, 16, 100)
            msg = "a"^n
            @test H8(msg) isa UInt8
            @test H16(msg) isa UInt16
        end
    end
end
