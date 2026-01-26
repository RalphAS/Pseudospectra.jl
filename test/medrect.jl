# medium rectangular (uses QR Lanczos)
using Pseudospectra, Test

@testset "Medium Rectangular" begin
    n = Pseudospectra.psathresholds.maxstdqr4hess + 10
    @show n
    G = Pseudospectra.grcar(n)
    # drop a column to test rectangular case
    A = G[:,1:end-1]
    # for a rectangular mtx, must provide axes
    opts = Dict{Symbol,Any}(:ax => [-1.0,3.0,-3.0,3.0],:npts=>60)

    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    @test ps_data.ps_dict[:algo] == :HessQR
    # non-Hessenberg case triggers an extra factorization
    A = G[:,1:end-1]
    A[n,1] = 1e-3
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    @test ps_data.ps_dict[:algo] == :HessQR
end
