# 60x59 med rectangular (uses QR Lanczos)
using Pseudospectra, Test

@testset "Medium Rectangular" begin
    G = Pseudospectra.grcar(60)
    # drop a column to test rectangular case
    A = G[:,1:end-1]
    # for a rectangular mtx, must provide axes
    opts = Dict{Symbol,Any}(:ax => [-1.0,3.0,-3.0,3.0],:npts=>60)

    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
end
