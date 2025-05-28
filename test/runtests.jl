using Test
using TinyCrypto

@testset "Curve tests" begin
    include("test_field.jl")
    include("test_montgomery.jl")
    include("test_edwards.jl")
    include("test_weierstrass.jl")
end

@testset "Example: DKG" begin
    include("../examples/dkg.jl")
    result = ExampleDKG.run()
    @test result == true
end
