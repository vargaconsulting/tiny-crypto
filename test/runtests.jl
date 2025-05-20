using Test
using TinyCrypto

@testset "Curve tests" begin
    @test 1 + 1 == 2
end

@testset "Example: DKG" begin
    include("../examples/dkg.jl")
    result = ExampleDKG.run()
    @test result == true
end
