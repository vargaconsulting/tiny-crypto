module Utils
export is_prime, primes, is_quadratic_residue, mod_inverse

    "Check if an integer is prime."
    function is_prime(n::T)::Bool where {T<:Integer}
        if n < 2
            return false
        end
        max = convert(T, floor(Int, sqrt(Int(n))))  # promote to Int for sqrt
        for i in 2:max
            if n % i == 0
                return false
            end
        end
        return true
    end

    "Return a vector of primes in the given integer range."
    function primes(range::UnitRange{T}) where {T<:Integer}
        return [n for n in range if is_prime(n)]
    end

    "Check if `a` is a quadratic residue mod `p` (assumes `p` is prime)."
    function is_quadratic_residue(a::T, p::T)::Bool where {T<:Integer}
        return a % p == 0 || powermod(a, (p - 1) ÷ 2, p) == one(T)
    end

    "Compute modular inverse using Fermat's Little Theorem (assumes prime `p`)."
    function mod_inverse(a::Integer, p::Integer)
        T = promote_type(typeof(a), typeof(p))
        a, p = T(a), T(p)
        return powermod(a, p - 2, p)
    end
end
