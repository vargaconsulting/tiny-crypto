module Utils
export is_prime, primes, is_quadratic_residue, mod_inverse, cofactor, prime_factors, is_square

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
    function max_prime_factor(N::Integer)
        n = N
        max_factor = 1
        d = 2
        while d * d <= n
            while n % d == 0
                max_factor = d
                n = div(n, d)
            end
            d += 1
        end
        return n > 1 ? n : max_factor
    end
    function cofactor(N::Integer)
        return div(N, max_prime_factor(N))
    end
    function prime_factors(N::Integer)::Vector{Int}
        factors = Int[]
        n = N
        d = 2
        while d * d <= n
            while n % d == 0
                push!(factors, d)
                n = div(n, d)
            end
            d += 1
        end
        if n > 1
            push!(factors, n)
        end
        return factors
    end
    function is_square(a::Integer, π::Integer)::Bool
        a % π == 0 && return true
        return powermod(a % π, (π - 1) >> 1, π) == 1
    end
end