module Field
export Fp, 𝔽ₚ
using Primes: isprime

    # Safe modular addition (no promotion, no overflow)
    function addmod(a::T, b::T, p::T) where {T<:Unsigned}
        sum = a + b
        if sum ≥ p || sum < a  # overflow or reduction needed
            sum -= p
        end
        return sum
    end

    # Safe modular multiplication using bitwise schoolbook method
    function mulmod(a::T, b::T, p::T) where {T<:Unsigned}
        result = zero(T)
        x, y = a, b  # <-- Use x instead of mutating a directly
        while y != 0
            if y & one(T) != 0
                result = addmod(result, x, p)
            end
            x = addmod(x, x, p)
            y >>= 1
        end
        return result
    end

    # Modular inverse using extended Euclidean algorithm
    function invmod(a::T, p::T) where {T<:Unsigned}
        a == 0 && throw(DivideError())
        t, newt = Int(0), Int(1)
        r, newr = Int(p), Int(a)

        while newr != 0
            q = r ÷ newr
            t, newt = newt, t - q * newt
            r, newr = newr, r - q * newr
        end

        t < 0 && (t += Int(p))
        return T(t)
    end

    struct Fp{T<:Unsigned, P}
        val::T
        function Fp{T,P}(x::Integer) where {T<:Unsigned, P}
            isprime(P) || throw(ArgumentError("Field modulus $P is not prime"))
            return new{T,P}(T(mod(x, P)))
        end
    end
    const 𝔽ₚ = Fp
    Base.Int64(x::Fp{T,P}) where {T<:Unsigned, P} = Int64(x.val)
    Base.Int32(x::Fp{T,P}) where {T<:Unsigned, P} = Int32(x.val)
    Base.UInt64(x::Fp{T,P}) where {T<:Unsigned, P} = UInt64(x.val)
    Base.UInt32(x::Fp{T,P}) where {T<:Unsigned, P} = UInt32(x.val)
    Base.Float64(x::Fp{T,P}) where {T<:Unsigned, P} = Float64(x.val)
    
    # Compact unicode display when P < 1000: e.g., 12𝔽₃₁, else fallback to 12𝔽ₚ
    function Base.show(io::IO, x::Fp{T,P}) where {T,P}
        val_str = "\e[1m$(x.val)\e[22m"  # ANSI bold
        if P < 1000000
            subscript = join(Char(0x2080 + d) for d in reverse(Base.digits(P)))
            print(io, "$(val_str)𝔽$subscript")
        else
            print(io, "$(val_str)⋅𝔽ₚ") 
        end
    end
    Base.show(io::IO, ::MIME"text/plain", x::Fp{T,P}) where {T,P} = show(io, x)

    # Promotion for interop
    Base.convert(::Type{Fp{T,P}}, x::Integer) where {T,P} = Fp{T,P}(x)
    Base.convert(::Type{Integer}, x::Fp{T,P}) where {T<:Unsigned, P} = x.val
    Base.convert(::Type{Int64}, x::Fp{T,P}) where {T<:Unsigned, P} = Int64(x.val)
    Base.convert(::Type{UInt64}, x::Fp{T,P}) where {T<:Unsigned, P} = UInt64(x.val)
    Base.convert(::Type{Int32}, x::Fp{T,P}) where {T<:Unsigned, P} = Int32(x.val)
    Base.convert(::Type{Int16}, x::Fp{T,P}) where {T<:Unsigned, P} = Int16(x.val)
    
    Base.promote_rule(::Type{Fp{T,P}}, ::Type{Int}) where {T<:Unsigned, P} = Fp{T,P}
    Base.promote_rule(::Type{Int}, ::Type{Fp{T,P}}) where {T<:Unsigned, P} = Fp{T,P}
    Base.promote_rule(::Type{Fp{T,P}}, ::Type{Integer}) where {T,P} = Fp{T,P}
    Base.promote_rule(::Type{Fp{T,P}}, ::Type{Fp{T,P}}) where {T,P} = Fp{T,P}

    # Basic arithmetic
    Base.mod(x::Fp{T,P}, y::Integer) where {T,P} = mod(x.val, y)
    Base.abs(x::Fp{T,P}) where {T,P} = x
    Base.inv(a::Fp{T,P}) where {T,P} = Fp{T,P}(invmod(a.val, T(P)))
    Base.:+(a::Fp{T,P}, b::Fp{T,P}) where {T,P} = Fp{T,P}(addmod(a.val, b.val, T(P)))
    Base.:+(a::Fp{T,P}, b::Integer) where {T,P} = a + Fp{T,P}(b)
    Base.:+(a::Integer, b::Fp{T,P}) where {T,P} = Fp{T,P}(a) + b
    Base.:-(a::Fp{T,P}) where {T,P} = Fp{T,P}(addmod(T(P) - a.val, zero(T), T(P)))
    Base.:-(a::Fp{T,P}, b::Fp{T,P}) where {T,P} = Fp{T,P}(addmod(a.val, T(P) - b.val, T(P)))
    Base.:-(a::Fp{T,P}, b::Integer) where {T,P} = a - Fp{T,P}(b)
    Base.:*(a::Fp{T,P}, b::Fp{T,P}) where {T,P} = Fp{T,P}(mulmod(a.val, b.val, T(P)))
    Base.:/(a::Fp{T,P}, b::Fp{T,P}) where {T,P} = a * inv(b)
    Base.:<<(x::Fp{T,P}, n::Integer) where {T,P} = Fp{T,P}(mod(x.val << n, T(P)))
    Base.:>>(x::Fp{T,P}, n::Integer) where {T,P} = Fp{T,P}(mod(x.val >> n, T(P)))
    Base.:-(a::Integer, b::Fp{T,P}) where {T,P} = Fp{T,P}(a) - b
    Base.:*(a::Integer, b::Fp{T,P}) where {T,P} = Fp{T,P}(a) * b
    Base.:*(a::Fp{T,P}, b::Integer) where {T,P} = a * Fp{T,P}(b)
    Base.:/(a::Integer, b::Fp{T,P}) where {T,P} = Fp{T,P}(a) / b
    Base.:/(a::Fp{T,P}, b::Integer) where {T,P} = a / Fp{T,P}(b)
    Base.:^(a::Fp{T,P}, n::Integer) where {T,P} = begin
        r = one(Fp{T,P})
        base = a
        while n > 0
            if n & 1 != 0
                r *= base
            end
            base *= base
            n >>= 1
        end
        return r
    end

    Base.zero(::Type{Fp{T,P}}) where {T,P} = Fp{T,P}(0)
    Base.one(::Type{Fp{T,P}}) where {T,P} = Fp{T,P}(1)
    Base.:(==)(a::Fp{T,P}, b::Fp{T,P}) where {T,P} = a.val == b.val
    Base.iszero(x::Fp{T,P}) where {T,P} = x.val == 0
    Base.isone(x::Fp{T,P}) where {T,P} = x.val == 1
    Base.getindex(x::Fp) = x.val
    Base.hash(x::Fp{T,P}, h::UInt) where {T,P} = hash((T, P, x.val), h)
    Base.iseven(x::Fp{T,P}) where {T,P} = iseven(x.val)
    Base.isodd(x::Fp{T,P}) where {T,P} = isodd(x.val)
    Base.copy(x::Fp{T,P}) where {T,P} = Fp{T,P}(x.val)
end # module
