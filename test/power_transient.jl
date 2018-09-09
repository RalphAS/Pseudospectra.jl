using Pseudospectra, Test

# Toeplitz matrix from S&P, Ch. 16
# Need a large domain for lower bound computation.
@testset "Matrix powers" begin
    opts = Dict{Symbol,Any}(:npts => 40,:ax => [-110,110,-110,110])
    A = diagm(1 => 100.0*ones(4)) + diagm(-2 => 1e-4*ones(3))
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)

    mtxpowersplot(gs,ps_data,20,lbmethod=:combined)
end
