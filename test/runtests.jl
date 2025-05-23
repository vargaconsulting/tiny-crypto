using Test
using TinyCrypto

@testset "Curve tests" begin
    include("test_hash.jl")
    include("test_field.jl")
end

@testset "Example: DKG" begin
    include("../examples/dkg.jl")
    result = ExampleDKG.run()
    @test result == true
end
