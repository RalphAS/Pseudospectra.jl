# To run only selected tests:
# $julia -e 'include("test/runtests.jl")' TEST1...

using Pseudospectra, Base.Test

tests = ["basic","arpack1","medcplx","medrect","smrect","spdirect",
          "projection_and_transient","radius_abscissa","negproj","eigsnoconv",
          "numrange", "linmap", "big"
          ]

if length(ARGS) > 0
    tests = ARGS
end

include("testplotter.jl")

for f in tests
    println("Testing $f functionality")
    fnam = f*".jl"
    include(fnam)
end
