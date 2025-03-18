#   Copyright © <2010-2025> Varga Consulting, Toronto, On     info@vargaconsulting.ca
macro attach(st, fields...)
	block = Expr(:block)
	for f in fields
		push!(block.args, :($f = $st.$f) )
    end
	esc(:($block))
end
macro define(T, fields...)
	block = Expr(:block)
	for f in fields
		push!(block.args, :($f ::$T ))
    end
	return esc(:($block))
end

function is_prime(n::Int)::Bool
    if n < 2 return false end
    for i in 2:floor(Int, sqrt(n))
        if n % i == 0
            return false
        end
    end
    return true
end

function primes(range::UnitRange{Int})
    return [n for n in range if is_prime(n)]
end

function is_quadratic_residue(a::Int, p::Int)::Bool
    if a % p == 0
        return true  # 0 is trivially a quadratic residue
    end
    result = powermod(a, (p - 1) ÷ 2, p)
    return result == 1
end
function mod_inverse(a::T, p::T) where T # Modular inverse
    return powermod(a, p - 2, p)         # Fermat's little theorem
end
