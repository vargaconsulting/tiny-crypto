using Test
using TinyCrypto

@testset "Curve tests" begin
    include("test_field.jl")
    include("test_montgomery.jl")
end

@testset "Example: DKG" begin
    include("../examples/dkg.jl")
    result = ExampleDKG.run()
    @test result == true
end
