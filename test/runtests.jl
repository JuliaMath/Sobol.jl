using Sobol
using Base.Test

# test integrand from Joe and Kuo paper (integrates to 1)
function testfunc(x::Vector{Float64})
    f = 1.0
    for j = 1:length(x)
        cj = j^0.3333333333333333333
        f *= (abs(4*x[j] - 2) + cj) / (1 + cj)
    end
    return f
end

function testint(N, n)
    s = SobolSeq(N)
    x = Array(Float64, N)
    skip(s, n)
    sum = 0.0
    for j = 1:n
        sum += testfunc(next!(s, x))
    end
    sum / n
end

# test results from Joe and Kuo (2003)
N = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1111]
n = [1009, 1997, 4001, 8009, 16001, 32003, 64007, 128021]
JoeKuo = [0.9743 0.9544 0.9838 0.9982 0.9995 0.9976 0.9969 0.9975;
          1.0157 0.9288 0.9649 0.9866 0.9985 0.9817 0.9955 0.9995;
          1.1109 0.8877 0.9070 1.0312 1.0216 0.9772 0.9744 0.9925;
          0.9770 0.8715 0.9305 0.9546 1.0060 0.9712 0.9770 0.9843;
          1.1780 0.8908 0.9710 0.9527 1.0203 0.9754 0.9779 0.9725;
          1.0567 0.8548 0.9700 0.9704 0.9768 0.9634 0.9745 0.9670;
          0.8735 0.7674 1.0130 0.9874 0.9677 0.9416 0.9580 0.9564;
          0.8752 0.8241 1.0067 1.0078 0.9987 0.9236 0.9517 0.9475;
          0.9139 0.8822 1.0926 1.0560 1.0210 0.9454 0.9330 0.9359;
          0.7779 1.0225 1.0593 1.0991 1.0470 0.9115 0.9688 0.9615;
          0.8348 1.0443 1.0710 1.0862 1.0483 0.9188 0.9735 0.9683;
          0.7524 0.9139 0.9592 0.9703 1.0116 0.9124 1.0411 0.9904]

for i = 1:length(N)
    println("testing integrals of dimension $(N[i])")
    for j = 1:length(n)
        println("  ... with $(n[j]) samples")
        t = testint(N[i], n[j])
        @test_approx_eq_eps t JoeKuo[i,j] 1e-4
    end
end
