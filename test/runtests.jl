using Sobol, Test

# compare results with results from C++ code sobol.cc published on
# http://web.maths.unsw.edu.au/~fkuo/sobol/
# with new-joe-kuo-6.21201 file used as input.
#
# Command line used to generate output is
#  ./sobol $N 1024 new-joe-kuo-6.21201 > exp_results_$N
#
# For each of the dimensions below
# (except 21201 where only 64 samples are generated)

@testset "published results" begin
    dimensions = [5, 10, 20, 50, 100, 500, 21201]

    for dim in dimensions
        println("Testing dimension $(dim)")
        open(joinpath(dirname(@__FILE__), "results", "exp_results_$(dim)")) do exp_results_file
            s = SobolSeq(dim)
            x = zeros(dim)
            for line in eachline(exp_results_file)
                values = [parse(Float64, item) for item in split(line)]
                if length(values) > 0
                    @test x == values
                    next!(s, x)
                end
            end
        end
    end
end

using Base.Iterators: take
@testset "iterators" begin
    # issue #8
    @test [x[1] for x in collect(take(Sobol.SobolSeq(1),5))] == [0.5,0.75,0.25,0.375,0.875]
end

@testset "scaled" begin
    # ScaledSobolSeq constructors
    lb = [-1,0,0]
    ub = [1,3,2]
    N = length(lb)
    s = SobolSeq(lb,ub)
    @test s isa ScaledSobolSeq{3,Int}
    @test eltype(s) == Vector{Float64}
    @test eltype(SobolSeq(Float32.(lb),Float32.(ub))) == Vector{Float32}
    @test first(s) == [0,1.5,1]
    @test first(SobolSeq((x for x in lb), (x for x in ub))) == [0,1.5,1]
    @test SobolSeq(N,lb,ub) isa ScaledSobolSeq{3,Int}
    @test_throws BoundsError SobolSeq(2,lb,ub)
end

@testset "skip exactly n" begin
    # issue #22
    # want to get the nth element
    n = 10

    s1 = SobolSeq(1)
    x1 = [next!(s1)[1] for i in 1:n][n]

    s2 = SobolSeq(1)
    skip(s2, n-1, exact=true)
    x2 = next!(s2)[1]

    @test x1 == x2
end
