module ECDSA
    import Base: sign    
    using ..TinyCrypto: H₈, invmod, isodd
    using ..TinyCrypto: Weierstrass, Point, ECPoint

    export genkey, sign, verify

    function genkey(curve::Weierstrass)
        priv = rand(1:curve.order - 1)
        pub  = priv * curve.G
        return priv, pub
    end

    function sign(curve::Weierstrass, priv, msg::AbstractString)
        z = Int(H₈(msg)) % curve.order  # hash message to scalar mod curve order

        while true
            k = rand(1:curve.order - 1)
            R = k * curve.G
            r = mod(R.point.x, curve.order)
            if r == 0 continue end  # Ensure r ≠ 0
            k_inv = invmod(k, curve.order)
            s = mod(k_inv * (z + r * priv), curve.order)
            if s ≠ 0
                parity = isodd(R.point.y.val) ? 1 : 0
                return (r = r, s = s, v = parity)
            end
        end
    end

    function verify(curve::Weierstrass, pub, sig, msg::AbstractString)
        r, s = sig.r, sig.s
        if !(1 <= r < curve.order && 1 <= s < curve.order) return false end
        #if !is_point_on_curve(pub, curve) return false end
        z = Int(H₈(msg)) % curve.order
        s_inv = invmod(s, curve.order)
        u₁ = z * s_inv % curve.order
        u₂ = r * s_inv % curve.order
        P = u₁ * curve.G + u₂ * pub
        return mod(P.point.x, curve.order) == r
    end
    function ecrecover(curve::Weierstrass, msg::AbstractString, r::Int, s::Int, v::Int)
        n, π = curve.order, curve.π
        z = Int(H₈(msg)) % n
        x_candidates = [r + k*n for k in 0:div(π,n) if r + k*n < π]
    
        for x in x_candidates
            α = mod(x^3 + curve.a * x + curve.b, π)
            y_candidates = [y for y in 0:π-1 if y^2 % π == α]
            for y in y_candidates
                if isodd(y) == v
                    R = ECPoint(x, y, curve)
                    r_inv = invmod(r, n)
                    sR = s * R
                    zG = z * curve.G
                    neg_zG = ECPoint(zG.point.x, -zG.point.y, curve)  # handle point negation properly
                    pub = r_inv * (sR + neg_zG)
                    return pub
                end
            end
        end
        return nothing
    end
    function ecrecover(curve::Weierstrass, msg::AbstractString, signature::NTuple{3,Int})
        r,s,v = signature
        return ecrecover(curve, msg, r,s,v)
    end
    function ecrecover(curve::Weierstrass, msg::AbstractString, signature::NamedTuple{(:r, :s, :v), <:Tuple{<:Integer, <:Integer, <:Integer}})
        return ecrecover(curve, msg, signature.r, signature.s, signature.v)
    end
end