using Test
using TinyCrypto

#FIXME: fails becuase of type casts, replace this with field arithmatics
@testset "ECPoint operator overloads and infinity behavior" begin
    C = Weierstrass(30:103, 0:2, 6:100)
    G = ECPoint(C.G, C)
    O = Infinity(C)

    # identity element
    @test O + G == G
    @test G + O == G
    @test 0 * G == O
    @test 1 * G == G
    @test 2 * G == G + G
    @test 3 * G == 2 * G + G

    # show output includes ∈
    io = IOBuffer()
    show(io, G)
    str = String(take!(io))
    @test occursin("∈", str)

    # show output for ∞
    #io_∞ = IOBuffer()
    #show(io_∞, infinity(T))
    #str_∞ = String(take!(io_∞))
    #@test str_∞ == "(∞,∞)"

    # is_infinity predicate
    @test is_infinity(O.point)
    @test !is_infinity(G.point)

    minus_G = -G
    @test G + minus_G == O
    @test 2 * G - G == G
    @test 2 * G == G + G              # point doubling
    @test G + 2 * G == 2 * G + G      # commutitativity
    @test (G + G) + G == G + (G + G)  # associativity
    @test -O == O                     # inverse
    @test O + O == O                  # 
    @test 2 * G + G == 3 * G          # mixed additions
    @test G + G + G == 3 * G          # ditto

    @test G == G                      # reflexivity/symmetry
    @test G + O == O + G

    @test 5 * G == 2 * G + 3 * G      # higher order scalar multiply

    @test typeof(2 * G) == typeof(G)
end
#=
@testset "Group law corner cases" begin
    C = Weierstrass(17:19, 1:3, 1:3)
    G = ECPoint(C.G, C)
    O = ECPoint(infinity(Int), C)

    @test 2 * G == G + G
    @test G + 2 * G == 2 * G + G
    @test (G + G) + G == G + (G + G)
    @test -O == O
    @test O + O == O
    @test 2 * G + G == 3 * G
    @test G + G + G == 3 * G
    @test -2 * G == -G + -G
end
=#