# verify that a LinearMap makes it through the chain

using Pseudospectra, LinearMaps, Base.Test

@testset "Linear map" begin

    A = LinearMap(Pseudospectra.convdiff_fd(10))

    # d,nconv,niter,nmult,resid = Pseudospectra.xeigs(A,I; nev=20, ncv=60,
    #                                                 which=:SR )

    opts = Dict{Symbol,Any}(:ARPACK_plots => false,:npts=>20)
    opts[:arpack_opts] = ArpackOptions{eltype(A)}(which=:SR,nev=6,ncv=60)

    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
end
