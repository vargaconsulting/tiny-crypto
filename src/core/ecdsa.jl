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
end