# complex dense, projected w/ transient
# Compare to figures in S&P, chap. 22.
# ∃ small discrepancies for large ϵ,
# but we are consistent w/ EigTool.
using Pseudospectra, Test

@testset "Projection/Transients" begin
    # for a sensible projection, must provide axes
    opts = Dict{Symbol,Any}(:ax => [-1,0.3,-1.2,0.2],:npts=>100,:proj_lev=>1.5)
    A = Pseudospectra.orrsommerfeld(200,10000,1.02)

    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)

    mtxexpsplot(gs,ps_data,1.0,100,lbmethod=:old)

    # additional plots to compare w/ book Fig. 22.6:

    # modeplot(ps_data,gs,0,0.0-0.2im)
    # modeplot(ps_data,gs,1,-0.1-0.8im)
end
