# test simple plotting interface

using Pseudospectra, Test

@testset "Undriven" begin
    A = Pseudospectra.grcar(40)
    spectralportrait(A)
end
