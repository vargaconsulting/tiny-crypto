module Hash
    export H, H₈,H₁₆, H8, H16

    const INIT8 = UInt8[0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0]
    const K8 = UInt8[0x1f, 0x2b, 0x3c, 0x4d, 0x5e, 0x6a, 0x7d, 0x8e]

    rotate(x::UInt8, n::Int) = ((x >> n) | (x << (8 - n))) % UInt8
    rotate(x::UInt16, n::Int) = ((x >> n) | (x << (16 - n))) % UInt16

    function pad8(msg::Vector{UInt8})
        len = length(msg)
        padded = copy(msg)
        push!(padded, 0x80)
        while (length(padded) % 8) != 7
            push!(padded, 0x00)
        end
        push!(padded, UInt8(len % 256))
        return padded
    end

    function compress_block8(block::Vector{UInt8}, state::Vector{UInt8})
        W = copy(block)
        for i in 1:8
            W[i] ⊻= K8[i]
            W[i] = rotate(W[i], i % 8 + 1)
        end
        for i in 1:8
            state[i] = rotate(state[i] ⊻ W[i], i)
        end
    end

    function tiny_hash8(msg::Vector{UInt8})
        padded = pad8(msg)
        state = copy(INIT8)
        for i in 1:8:length(padded)
            block = padded[i:i+7]
            compress_block8(block, state)
        end
        result = state[1]
        for i in 2:8
            result ⊻= rotate(state[i], i)
        end
        return result
    end
    tiny_hash8(msg::String) = tiny_hash8(Vector{UInt8}(codeunits(msg)))

    function tiny_hash16(msg::Vector{UInt8})
        acc = UInt16(0xabcd)
        for (i, byte) in enumerate(msg)
            acc ⊻= UInt16(byte) << (i % 8)
            acc = (acc << 5) ⊻ (acc >> 3)
        end
        return acc
    end
    tiny_hash16(msg::String) = tiny_hash16(Vector{UInt8}(codeunits(msg)))

    function H8(msg::Union{String, Vector{UInt8}})
        bytes = msg isa String ? Vector{UInt8}(codeunits(msg)) : msg
        return tiny_hash8(bytes)
    end

    function H16(msg::Union{String, Vector{UInt8}})
        bytes = msg isa String ? Vector{UInt8}(codeunits(msg)) : msg
        return tiny_hash16(bytes)
    end

    H(msg::Union{String, Vector{UInt8}}) = H8(msg)
    H₈(msg::Union{String, Vector{UInt8}}) = H8(msg)
    H₁₆(msg::Union{String, Vector{UInt8}}) = H16(msg)
end