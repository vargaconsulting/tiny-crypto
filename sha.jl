module SHA8
export sha8

    const INIT = UInt8[0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0]
    const K = [0x1f, 0x2b, 0x3c, 0x4d, 0x5e, 0x6a, 0x7d, 0x8e]

    rotate(x::UInt8, n::Int) = ((x >> n) | (x << (8 - n))) % UInt8

    function pad(msg::Vector{UInt8})
        len = length(msg)
        padded = copy(msg)
        push!(padded, 0x80)  # add 1 bit

        while (length(padded) % 8) != 7
            push!(padded, 0x00)
        end

        push!(padded, UInt8(len % 256))  # length mod 256
        return padded
    end

    function compress_block(block::Vector{UInt8}, state::Vector{UInt8})
        W = copy(block)
        for i in 1:8
            W[i] ⊻= K[i]
            W[i] = rotate(W[i], i % 8 + 1)
        end

        for i in 1:8
            state[i] = rotate(state[i] ⊻ W[i], i)
        end
    end

    function sha8(msg::Vector{UInt8})
        padded = pad(msg)
        state = copy(INIT)

        for i in 1:8:length(padded)
            block = padded[i:i+7]
            compress_block(block, state)
        end

        # Fold the 8-byte state to a single byte
        result = state[1]
        for i in 2:8
            result ⊻= rotate(state[i], i)
        end

        return result
    end
    function sha8(msg::String)
        return sha8(collect(UInt8, msg))
    end

end # module

