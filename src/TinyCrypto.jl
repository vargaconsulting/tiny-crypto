module TinyCrypto

include("macroes.jl")
include("utils.jl")
include("core/ecpoint.jl")
include("curves/weierstrass.jl")

    using .Macroes
    using .Utils

    # Re-export core symbols for public API
    export @define, @attach
    export is_prime, primes, mod_inverse
    export Point, ∞, Curve, ECPoint
    export Weierstrass
    export infinity, is_infinity
    export point_add, scalar_mult, curve_points, is_generator, is_point_on_curve
end