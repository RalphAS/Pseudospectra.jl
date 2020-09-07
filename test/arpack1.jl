using Pseudospectra, Test

@testset "Sparse-ARPACK" begin
    # 400x400 real sparse
    A = Pseudospectra.convdiff_fd(10)
    #opts = Dict{Symbol,Any}(:ARPACK_plots => false,:npts=>40)
    opts = Dict{Symbol,Any}(:ARPACK_plots => true,:npts=>144,
                            :ax=>[30,90,-70,70])
    opts[:arpack_opts] = ArpackOptions{eltype(A)}(which=:SR,nev=20,ncv=60)

    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
end
