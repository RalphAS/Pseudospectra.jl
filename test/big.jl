# test on generic type

using Pseudospectra, Test, GenericSVD

@testset "Generic(Big)" begin

    @info("This test is expected to comment about lack of methods for BigFloat eigvals")
    A = Matrix{BigFloat}(Pseudospectra.grcar(8))
    # ax is needed until ∃ eigvals(BigFloat)
    opts = Dict{Symbol,Any}(:ax => [-1,3,-3,3], :npts => 20)
    ps_data = (@test_logs (:warn, r"^Failed to compute eigenvalues") new_matrix(A,opts))
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)

    # Just big enough to get to the inverse-Lanczos branch
    # this is a stunt, so don't waste time with usual npts.
    A = Matrix{BigFloat}(Pseudospectra.grcar(56))
    opts = Dict{Symbol,Any}(:ax => [-1,3,-3,3], :npts => 10)
    @test_logs (:warn, r"^Failed to compute eigenvalues") ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
end
