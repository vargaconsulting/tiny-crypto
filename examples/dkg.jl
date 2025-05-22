module ExampleDKG
using TinyCrypto, LinearAlgebra, Random

    function mod_inverse(a, p)
        return powermod(a, p - 2, p)
    end

    function lagrange_basis(x, j, 𝔽ₚ)
        xⱼ, x₀, N, D = x[j], 0, 1, 1
        for (k, xₖ) ∈ enumerate(x)
            if k != j
                N = (N * (x₀ - xₖ)) % 𝔽ₚ
                D = (D * (xⱼ - xₖ)) % 𝔽ₚ
            end
        end
        return (N * mod_inverse(D, 𝔽ₚ)) % 𝔽ₚ
    end

    # Compute the secret using Lagrange interpolation
    function y_intercept(x, shards, 𝔽ₚ)
        secret = 0
        for (j, s) ∈ enumerate(shards)
            ℓⱼ = lagrange_basis(x, j, 𝔽ₚ)
            secret = mod(secret + s * ℓⱼ, 𝔽ₚ)
        end
        return secret
    end

    function run()
        # Distributed Secret Construction Protocol with t treshold out of p participants  
        # Participants p:=7 (P₁,P₂,P₃,..,P₇) Threshold: t:=3 => cubic polinomial,  Finite Field Prime: 97 
        # Collaboratively construct a secret without any participant holding the full secret.
        p,t, 𝔽ₚ = 7, 3, 97 # p := participants, t:=threshold, 𝔽ₚ = field prime 

        # Step 1: Participants Generate Random Polynomials in the form: Qᵢ(x) = rᵢ + aᵢx + aᵢx² + .. + aᵢxᵗ⁻¹ 
        Q = rand(1:𝔽ₚ - 1, p, t + 1) #  Q ~ Uniform(1, 2²⁵⁶) a cryptographic uniform distribution
        # upto 𝔽ₚ - field prime for brevity we are using a small filed prime 97
        # such that rᵢ = Q[:,1] which is the first column, and the coeffs are: a = Q[:,2:end]
        Pₖ = mod(sum(Q[:,1]), 97)  # the sum of y intercepts is our private key constructed in a distributed fashion.

        # Step 2: Evaluate the random polynomials at same number of points as participants, so we can redistribute the  results ALL to ALL
        polynomial(degree) = x -> [x^i for i in 0:degree]
        f = polynomial(t)  # precompute the indeterminate values with degrees {0,1,2}
        # evaluate the indeterminates, the polynomial terms without the coefficients, at random points:
        x = randperm(𝔽ₚ - 1)[1:p]     # out of a unique public  set of x locations, distributed and used for all participants
        X = mod.(hcat(f.(x)...)', 𝔽ₚ) # we are to evalue the intederminates of the polynomial forming X Vandermonde matrix over 𝔽ₚ
        # then multiply the indetermintes with the Q coefficients, then  sum them up, obtaining Mᵉ[p x p] evaluated matrix
        # before sharing, each row belong to a single Pᵢ participant: 
        Mᵉ = mod.(X * Q', 𝔽ₚ)

        # Step 3: Each participants distribute their computed Mᵉ evaluated sahreable  vector to all other participants so that:
        # M[1,2] => P[2,1], M[1,3] => P[3,1] ... M[1,n] => P[n,1] ALL to ALL from the pattern this is just the matrix transpose of Mᵉ : 
        Mᵈ = Mᵉ'
        # Step 4: compute the shards, which is the rowsome of the received values on our 𝔽ₚ prime field plus the first element
        # each Pᵢ kepts for itself:
        S = mod.(sum(Mᵈ, dims=1), 𝔽ₚ)'
        # At this point we have S[i] shards constructed in a fashion that no single Pᵢ participants could possible reconstruct the 
        # random value, which happenes to be the PRIVATE KEY (a cryptographic number over a 𝔽ₚ prime field ), until this point
        # no discrete logarithm was involved, therefore Shors algorithm running on ( future) Quantum Computers doesn't ruin our day. 


        ## PROOF of RECONSTRUCTION  (should not be performed in production system, instead see homomorphic operations for signing ) ##
        # we are to reconstruct the y intercept of Q polynomial at x = 0, (y axis) from t + 1 'S' shards using Lagrange Interpolation:
        # s  = P(0) = Σ Sⱼ​⋅ℓⱼ(0) mod 𝔽ₚ where ℓⱼ(0) is the Langrange basis evaluated at 0 (zero) 
        # Luckily this is equivalent to y = V'c where y := shards, and V is the Vandermonde matrix constructed from the public
        # evaluation points, and we want c coefficients: 

        # V⁻¹ = adj(V) × (det V)⁻¹ mod 𝔽ₚ  (however in practice the Vandermonde matrix can be ill conditioned: do not use)
        # from Fermat little theorem we have a⁻¹ ≡ aᵖ⁻²  in other words:
        mod_inverse(a, p) = powermod(a, p - 2, p)


        # for arbitrary combinations  of the shards and public x nodes the y intercept 
        # is invariant, and happens to be the constructed private key
        i = sort(randperm(p)[1:t+1])
        pₖ = y_intercept(x[i], S[i], 𝔽ₚ)
        pₖ == Pₖ 

        ## PROACTIVE REFRESH: adding a random polynomial with y intercept being 0, does not change the private key
        Qᵤ = hcat(zeros(Int, p),rand(1:𝔽ₚ - 1, p, t)) 
        Yᵤ = sum(Qᵤ[:,1])                  # the y intercept should be zero
        Uᵉ = mod.(X * Qᵤ', 𝔽ₚ)             # the polynomials evaluated at x
        Uᵈ = Uᵉ'                           # distribute the values
        Sᵤ = mod.(sum(Uᵈ, dims=1), 𝔽ₚ)'    # compute the shard updates
        Sₙ = S + Sᵤ                        # do the update 
        i = sort(randperm(p)[1:t+1])       # pick t threshold, or more random participants
        pᵤ = y_intercept(x[i], Sₙ[i], 𝔽ₚ)  # recover the secret (this should not be done, instead use homomorphic computation)
        pᵤ == Pₖ                           # this should be equal
    end
end