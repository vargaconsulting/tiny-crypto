module Macroes

export @define, @attach

    macro attach(st, fields...)
        block = Expr(:block)
        for f in fields
            push!(block.args, :($f = $st.$f))
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
end