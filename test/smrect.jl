# small rectangular
using Pseudospectra, Test

@testset "Small Rectangular" begin
    # Grcar matrix is Hessenberg
    n = Pseudospectra.psathresholds.maxstdqr4hess - 10
    G = Pseudospectra.grcar(n)
    A = G[:,1:end-1]
    # for a rectangular mtx, must provide axes
    opts = Dict{Symbol,Any}(:ax => [-1.0,3.0,-3.0,3.0],:npts=>60)

    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    @test ps_data.ps_dict[:algo] == :rect_qr
    # non-Hessenberg case triggers an extra factorization
    A = G[:,1:end-1]
    A[n,1] = 1e-3
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    @test ps_data.ps_dict[:algo] == :rect_qz
end
