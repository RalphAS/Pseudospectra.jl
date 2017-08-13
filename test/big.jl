# test on generic type

using Pseudospectra, Base.Test, GenericSVD

@testset "Generic(Big)" begin

    A = Matrix{BigFloat}(Pseudospectra.grcar(8))
    # ax is needed until ∃ eigvals(BigFloat)
    opts = Dict{Symbol,Any}(:ax => [-1,3,-3,3], :npts => 20)
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)

    # Just big enough to get to the inverse-Lanczos branch
    # this is a stunt, so don't waste time with usual npts.
    A = Matrix{BigFloat}(Pseudospectra.grcar(56))
    opts = Dict{Symbol,Any}(:ax => [-1,3,-3,3], :npts => 10)
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
end
