#=

    Data on the primitive binary polynomials (a) and the corresponding
    starting values m, for Sobol sequences in up to 21201 dimensions,
    taken from:
        S. Joe and F. Y. Kuo, Constructing Sobol sequences with better two-dimensional projections,
        SIAM J. Sci. Comput. 30, 2635-2654 (2008).

    Figures used were derived from files provided on Joe Kuo's website
    http://web.maths.unsw.edu.au/~fkuo/sobol/ on 16/10/2010
    and are provided under the license below.

    Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.

        * Neither the names of the copyright holders nor the names of the
          University of New South Wales and the University of Waikato
          and its contributors may be used to endorse or promote products derived
          from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#

using DelimitedFiles

const sobol_a_file = joinpath(@__DIR__, "sobol_a.csv")
const sobol_minit_file = joinpath(@__DIR__, "sobol_minit.csv")

# we need to re-compile if these files change:
include_dependency(sobol_a_file)
include_dependency(sobol_minit_file)

#successive primitive binary-coefficient polynomials p(z)
#   = a_0 + a_1 z + a_2 z^2 + ... a_31 z^31, where a_i is the
#     i-th bit of sobol_a[j] for the j-th polynomial.
const sobol_a = vec(readdlm(sobol_a_file, ',', UInt32))

# starting direction #'s m[i] = sobol_minit[i][j] for i=0..d of the
# degree-d primitive polynomial sobol_a[j].
const sobol_minit = readdlm(sobol_minit_file, ',', UInt32)
