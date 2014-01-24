module Sobol

import NLopt.libnlopt
import Base: ndims, length, skip, next, start, done, show, writemime
export SobolSeq, ScaledSobolSeq, next!

type SobolSeq{N}
    p::Ptr{Void} # nlopt_sobol pointer
    function SobolSeq()
        (N < 1 || N > 1111) && error("invalid Sobol dimension")
        p = ccall((:nlopt_sobol_create,libnlopt), Ptr{Void}, (Cuint,), N)
        p == C_NULL && error("error initializing Sobol sequence")
        s = new(p)
        finalizer(s, destroy)
        s
    end
end
SobolSeq(N::Integer) = SobolSeq{N}()

destroy(s::SobolSeq) = 
    ccall((:nlopt_sobol_destroy,libnlopt), Void, (Ptr{Void},), s.p)

ndims{N}(s::SobolSeq{N}) = N

function show{N}(io::IO, s::SobolSeq{N})
    print(io, "$N-dimensional Sobol sequence on [0,1]^$N")
end
function writemime{N}(io::IO, ::MIME"text/html", s::SobolSeq{N})
    print(io, "$N-dimensional Sobol sequence on [0,1]<sup>$N</sup>")
end

function skip!{N}(s::SobolSeq{N}, n::Integer, x::Vector{Float64})
    (n < 0 || length(x) < N) && throw(BoundsError())
    ccall((:nlopt_sobol_skip,libnlopt), Void, (Ptr{Void}, Cuint, Ptr{Float64}),
          s.p, n, x)
    x
end

function skip{N}(s::SobolSeq{N}, n::Integer) 
    skip!(s, n, Array(Float64, N))
    nothing
end

function next!{N}(s::SobolSeq{N}, x::Vector{Float64})
    length(x) < N && throw(BoundsError())
    ccall((:nlopt_sobol_next01,libnlopt), Void,
          (Ptr{Void}, Ptr{Float64}), s.p, x)
    x
end
next{N}(s::SobolSeq{N}) = next!(s, Array(Float64,N))

function next!{N}(s::SobolSeq{N}, x::Vector{Float64},
                  lb::Vector, ub::Vector)
    length(x) < N && throw(BoundsError())
    ccall((:nlopt_sobol_next,libnlopt), Void,
          (Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
          s.p, x, float64(lb), float64(ub))
    x
end
next{N}(s::SobolSeq{N}, lb::Vector, ub::Vector) =
    next!(s, Array(Float64,N), lb, ub)

# Make an iterator so that we can do "for x in SobolSeq(...)".
# Technically, the Sobol sequence ends after 2^32-1 points, but NLopt
# falls back on pseudorandom numbers after this.  In practice, one is
# unlikely to reach that point.
start(s::SobolSeq) = s
next(s::SobolSeq, s_::SobolSeq) = (next(s), s_)
done(s::SobolSeq, s_::SobolSeq) = false

# Convenience wrapper for scaled Sobol sequences

type ScaledSobolSeq{N}
    s::SobolSeq{N}
    lb::Vector{Float64}
    ub::Vector{Float64}
    ScaledSobolSeq(lb::Vector{Float64}, ub::Vector{Float64}) =
        new(SobolSeq(N), lb, ub)
end
SobolSeq(N::Integer, lb, ub) =
    ScaledSobolSeq{N}(copy!(Array(Float64,N),lb), copy!(Array(Float64,N),ub))

next!(s::ScaledSobolSeq, x::Vector{Float64}) = next!(s.s, x, s.lb, s.ub)
next{N}(s::ScaledSobolSeq{N}) = next!(s, Array(Float64,N))

start(s::ScaledSobolSeq) = s
next(s::ScaledSobolSeq, s_::ScaledSobolSeq) = (next(s), s_)
done(s::ScaledSobolSeq, s_::ScaledSobolSeq) = false

skip(s::ScaledSobolSeq, n) = skip(s.s, n)

function show{N}(io::IO, s::ScaledSobolSeq{N})
    lb = s.lb; ub = s.ub
    print(io, "$N-dimensional scaled Sobol sequence on [$(lb[1]),$(ub[1])]")
    cnt = 1
    for i = 2:N
        if lb[i] == lb[i-1] && ub[i] == ub[i-1]
            cnt += 1
        else
            if cnt > 1
                print(io, "^", cnt)
            end
            print(io, " x [", lb[i], ",", ub[i], "]")
            cnt = 1
        end
    end
    if cnt > 1
        print(io, "^", cnt)
    end
end

function writemime{N}(io::IO, ::MIME"text/html", s::ScaledSobolSeq{N})
    lb = s.lb; ub = s.ub
    print(io, "$N-dimensional scaled Sobol sequence on [$(lb[1]),$(ub[1])]")
    cnt = 1
    for i = 2:N
        if lb[i] == lb[i-1] && ub[i] == ub[i-1]
            cnt += 1
        else
            if cnt > 1
                print(io, "<sup>", cnt, "</sup>")
            end
            print(io, " × [", lb[i], ",", ub[i], "]")
            cnt = 1
        end
    end
    if cnt > 1
        print(io, "<sup>", cnt, "</sup>")
    end
end

end # module
