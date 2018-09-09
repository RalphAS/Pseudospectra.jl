# test basic case, with mode plots

using Pseudospectra, Test

@testset "Basic" begin

    A = Pseudospectra.grcar(40)

    opts = Dict{Symbol,Any}()
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    modeplot(ps_data,gs,0,0.5+2.0im)
    modeplot(ps_data,gs,1,0.5+2.0im)

end
