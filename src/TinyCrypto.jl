module TinyCrypto

include("macroes.jl")
include("utils.jl")
include("field.jl")

using .Macroes
using .Field
using .Utils

include("core/ecpoint.jl")
include("curves/weierstrass.jl")
include("curves/montgomery.jl")

    # Re-export core symbols for public API
    export @define, @attach
    export Fp, 𝔽ₚ
    export is_prime, primes, mod_inverse
    export Point, ∞, Curve, ECPoint, Infinity
    export Weierstrass, Montgomery
    export infinity, is_infinity
    export point_add, scalar_mult, curve_points, subgroup_points, is_generator, is_point_on_curve, inverse, is_singular
end

