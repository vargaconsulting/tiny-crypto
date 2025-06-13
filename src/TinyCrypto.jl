module TinyCrypto

include("macroes.jl")
include("utils.jl")
include("field.jl")
include("core/hash.jl")

using .Macroes
using .Field
using .Utils
using .Hash

include("core/ecpoint.jl")
include("curves/weierstrass.jl")
include("curves/montgomery.jl")
include("curves/edwards.jl")
    # Re-export core symbols for public API
    export @define, @attach
    export H, H₈, H₁₆, H8, H16
    export Fp, 𝔽ₚ
    export is_prime, primes, mod_inverse
    export Point, ∞, Curve, AffinePoint, Infinity
    export Weierstrass, Montgomery, Edwards, TwistedEdwards
    export infinity, is_infinity, is_identity
    export point_add, scalar_mult, curve_points, subgroup_points, is_generator, is_point_on_curve, inverse, is_singular


    get_point(p::AffinePoint) = p.point
    get_point(p::Point) = p

    function print_point_collection(io::IO, items, prefix::AbstractString, postfix::AbstractString; max_display=3)
        get_coords(p) = is_infinity(p) ? "(∞,∞)" : "($(get_point(p).x.val),$(get_point(p).y.val))"
        # Support Dict-style formatting
        point_strs = if isa(items, AbstractDict)
            [string(get_coords(k), " ↦ ", v) for (k, v) in items]
        else
            [get_coords(p) for p in items]
        end

        max_width = displaysize(stdout)[2]
        total_width = 2 + sum(length, point_strs) + (length(point_strs) - 1) * 2

        display_points = if total_width > max_width && length(point_strs) > 2 * max_display
        vcat(point_strs[1:max_display], "…", point_strs[end - max_display + 1:end])
        else
        point_strs
        end

        print(io, prefix, " {", join(display_points, ", "), "} ⊂ ", postfix)
    end

    function Base.show(io::IO, P::Point{F}) where {F<:Fp}
        if is_infinity(P)
            print(io, "(∞,∞)")
        else
            print(io, "(", P.x, ",", P.y, ")")
        end
    end

    subscript(x::Integer) = join(Char(0x2080 + d) for d in reverse(digits(x)))
    Base.show(io::IO, ::Type{Fp{T,P}}) where {T, P} = print(io, "𝔽", subscript(P))
    Base.show(io::IO, x::Fp{T,P}) where {T,P} = print(io, "\e[1m$(x.val)\e[22m")
    Base.show(io::IO, ::MIME"text/plain", x::Fp{T,P}) where {T,P} = show(io, x)


    Base.show(io::IO, P::AffinePoint{T,C}) where {T,C} = print(io, "$(P.point) ∈ ", typeof(P.curve))
    Base.show(io::IO, ::MIME"text/plain", s::Vector{AffinePoint{F,C}}) where {F,C<:Curve} = print_point_collection(io, s, "vector of", string(C))
    Base.show(io::IO, ::MIME"text/plain", s::Set{AffinePoint{F,C}}) where {F,C<:Curve} = print_point_collection(io, s, "set of", string(C))
    Base.show(io::IO, ::MIME"text/plain", v::Vector{Point{F}}) where {F<:Fp} = print_point_collection(io, v, "vector of", string(F))
    Base.show(io::IO, ::MIME"text/plain", v::Set{Point{F}}) where {F<:Fp} = print_point_collection(io, v, "set of" , string(F))
    Base.show(io::IO, ::MIME"text/plain", d::Dict{AffinePoint{F,C},V}) where {F,C<:Curve, V} = print_point_collection(io, d, "dict of", string(C))
    Base.show(io::IO, ::MIME"text/plain", d::Dict{Point{F},V}) where {F<:Fp, V} = print_point_collection(io, d, "dict of", string(F))
    Base.show(io::IO, ::MIME"text/plain", t::Tuple{Vararg{AffinePoint{F,C}}}) where {F,C<:Curve} = print_point_collection(io, t, "tuple of", string(C))
    Base.show(io::IO, ::MIME"text/plain", t::Tuple{Vararg{Point{F}}}) where {F<:Fp} = print_point_collection(io, t, "tuple of", string(F))

    function Base.show(io::IO, v::Vector{AffinePoint{F,C}}) where {F,C<:Curve} 
        print(io, "AffinePoint{$(F), $(C)}[")
        print_point_collection(io, v, ""; max_display=3)
        print(io, "]")
    end
end
