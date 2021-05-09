# test consistency of multithreading option

using Pseudospectra, LinearAlgebra

@testset "Multithreading" begin
@info "Testing multihreading with nt=$(Threads.nthreads())"
for n in [20,60]
    A = randn(n,n) + 1im * randn(n,n)
    Tschur,U,eigA  = schur(A)
    ax = [-3,3,-3,3]
    npts = 100
    opts = Dict{Symbol,Any}()
    Z0,x,y,levels,err,Tproj,eigAproj = psa_compute(Tschur,npts,
                                                   ax,eigA,opts)
    opts = Dict{Symbol,Any}(:threaded => true)
    Z1,x,y,levels,err,Tproj,eigAproj = psa_compute(Tschur,npts,
                                                   ax,eigA,opts)
    @test norm(Z0 - Z1) / norm(Z0) < 1e-4
end
end
