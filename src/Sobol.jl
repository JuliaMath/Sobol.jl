VERSION < v"0.7.0-beta2.199" && __precompile__()

module Sobol
using Compat, Compat.Random
export SobolSeq, ScaledSobolSeq, next!, next
@static if isdefined(Base, :next)
    import Base: next
end

include("soboldata.jl") #loads `sobol_a` and `sobol_minit`

abstract type AbstractSobolSeq{N} end

# N iis the dimension of sequence being generated
mutable struct SobolSeq{N} <: AbstractSobolSeq{N}
    m::Array{UInt32,2} #array of size (sdim, 32)
    x::Array{UInt32,1} #previous x = x_n, array of length sdim
    b::Array{UInt32,1} #position of fixed point in x[i] is after bit b[i]
    n::UInt32 #number of x's generated so far
end

ndims(s::AbstractSobolSeq{N}) where {N} = N::Int

function SobolSeq(N::Int)
    (N < 0 || N > (length(sobol_a) + 1)) && error("invalid Sobol dimension")

    m = ones(UInt32, (N, 32))

    #special cases
    N == 0 && return(SobolSeq{0})(m,UInt32[],UInt32[],zero(UInt32))
    #special cases 1
    N == 1 && return(SobolSeq{N}(m,UInt32[0],UInt32[0],zero(UInt32)))

    for i = 2:N
        a = sobol_a[i-1]
        d = floor(Int, log2(a)) #degree of poly

        #set initial values of m from table
        m[i, 1:d] = sobol_minit[1:d, i - 1]
        #fill in remaining values using recurrence
        for j = (d+1):32
            ac = a
            m[i,j] = m[i,j-d]
            for k = 0:d-1
                @inbounds m[i,j] = m[i,j] ⊻ (((ac & one(UInt32)) * m[i, j-d+k]) << (d-k))
                ac >>= 1
            end
        end
    end
    SobolSeq{N}(m,zeros(UInt32,N),zeros(UInt32,N),zero(UInt32))
end
SobolSeq(N::Integer) = SobolSeq(Int(N))

# 1/2^m for m = 1:32
const scale2m = exp2.(-(1:32))

function next!(s::SobolSeq, x::AbstractVector{<:AbstractFloat})
    length(x) != ndims(s) && throw(BoundsError())

    if s.n == typemax(s.n)
        return rand!(x)
    end

    s.n += one(s.n)
    c = UInt32(trailing_zeros(s.n))
    sb = s.b
    sx = s.x
    sm = s.m
    for i=1:ndims(s)
        @inbounds b = sb[i]
        if b >= c
            @inbounds sx[i] = sx[i] ⊻ (sm[i,c+1] << (b-c))
            @inbounds x[i] = sx[i] * scale2m[b+1]
        else
            @inbounds sx[i] = (sx[i] << (c-b)) ⊻ sm[i,c+1]
            @inbounds sb[i] = c
            @inbounds x[i] = sx[i] * scale2m[c+1]
        end
    end
    return x
end
next!(s::SobolSeq) = next!(s, Array{Float64,1}(undef, ndims(s)))

@deprecate next(s::AbstractSobolSeq) next!(s)

# if we know in advance how many points (n) we want to compute, then
# adopt the suggestion of the Joe and Kuo paper, which in turn
# is taken from Acworth et al (1998), of skipping a number of
# points equal to the largest power of 2 smaller than n
function skip!(s::SobolSeq, n::Integer, x)
    nskip = 1 << floor(Int,log2(n))
    for unused=1:nskip; next!(s,x); end
    return nothing
end
Base.skip(s::SobolSeq, n::Integer) = skip!(s, n, Array{Float64,1}(undef, ndims(s)))

function Base.show(io::IO, s::SobolSeq)
    print(io, "$(ndims(s))-dimensional Sobol sequence on [0,1]^$(ndims(s))")
end
function Base.show(io::IO, ::MIME"text/html", s::SobolSeq)
    print(io, "$(ndims(s))-dimensional Sobol sequence on [0,1]<sup>$(ndims(s))</sup>")
end

# Make an iterator so that we can do "for x in SobolSeq(...)".
# Technically, the Sobol sequence ends after 2^32-1 points, but it
# falls back on pseudorandom numbers after this.  In practice, one is
# unlikely to reach that point.
@static if isdefined(Base, :iterate)
    Base.iterate(s::AbstractSobolSeq, state=nothing) = (next!(s), state)
else
    Base.start(s::AbstractSobolSeq) = nothing
    Base.next(s::AbstractSobolSeq, state) = (next!(s), state)
    Base.done(s::AbstractSobolSeq, state) = false
end
Base.eltype(::Type{<:AbstractSobolSeq}) = Vector{Float64}
Compat.IteratorSize(::Type{<:AbstractSobolSeq}) = Base.IsInfinite()
Compat.IteratorEltype(::Type{<:AbstractSobolSeq}) = Base.HasEltype()

# Convenience wrapper for scaled Sobol sequences

struct ScaledSobolSeq{N} <: AbstractSobolSeq{N}
    s::SobolSeq{N}
    lb::Vector{Float64}
    ub::Vector{Float64}
    function ScaledSobolSeq{N}(lb::Vector{<:Real}, ub::Vector{<:Real}) where {N}
        length(lb)==length(ub)==N || throw(DimensionMismatch("lb and ub do not have length $N"))
        new(SobolSeq(N), lb, ub)
    end
end
SobolSeq(N::Integer, lb, ub) =
    ScaledSobolSeq{Int(N)}(copyto!(Vector{Float64}(undef,N), lb), copyto!(Vector{Float64}(undef,N), ub))
SobolSeq(lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real}) =
    SobolSeq(length(lb), lb, ub)


function next!(s::SobolSeq, x::AbstractVector{<:AbstractFloat},
               lb::AbstractVector, ub::AbstractVector)
    length(x) < ndims(s) && throw(BoundsError())
    next!(s,x)
    for i=1:ndims(s)
        x[i] = lb[i] + (ub[i]-lb[i]) * x[i]
    end
    return x
end
next!(s::SobolSeq{N}, lb::AbstractVector, ub::AbstractVector) where {N} = next!(s, Vector{Float64}(undef, N), lb, ub)

next!(s::ScaledSobolSeq, x::AbstractVector{<:AbstractFloat}) = next!(s.s, x, s.lb, s.ub)
next!(s::ScaledSobolSeq) = next!(s.s, Array{Float64,1}(undef, ndims(s)), s.lb, s.ub)
next(s::ScaledSobolSeq) = next!(s, Vector{Float64}(undef, ndims(s)))

Base.skip(s::ScaledSobolSeq, n) = skip(s.s, n)

function Base.show(io::IO, s::ScaledSobolSeq)
    lb = s.lb; ub = s.ub
    N = ndims(s)
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

function Base.show(io::IO, ::MIME"text/html", s::ScaledSobolSeq)
    N = ndims(s)
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
