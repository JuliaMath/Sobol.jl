using Sobol, Compat
using Compat.Test

# compare results with results from C++ code sobol.cc published on
# http://web.maths.unsw.edu.au/~fkuo/sobol/
# with new-joe-kuo-6.21201 file used as input.
#
# Command line used to generate output is
#  ./sobol $N 1024 new-joe-kuo-6.21201 > exp_results_$N
#
# For each of the dimensions below
# (except 21201 where only 64 samples are generated)

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

# issue #8
using Base.Iterators: take
@test [x[1] for x in collect(take(Sobol.SobolSeq(1),5))] == [0.5,0.75,0.25,0.375,0.875]

# ScaledSobolSeq constructors
lb = [-1.0,0,0]
ub = [1.0,3,2]
N = length(lb)
@test typeof(ScaledSobolSeq(N,lb,ub)) == typeof(SobolSeq(N,lb,ub))
@test typeof(ScaledSobolSeq(N,lb,ub)) == typeof(ScaledSobolSeq(lb,ub))
