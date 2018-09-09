using Pseudospectra, Test

@testset "Zoom" begin

    A = Pseudospectra.grcar(40)
    opts = Dict{Symbol,Any}()
    ps_data = new_matrix(A,opts)

    driver!(ps_data,opts,gs)

    zoomin!(gs,ps_data,zkw=1.2+1.8im)
    @test ps_data.zoom_pos == 2

    # this should return to original portrait
    zoomout!(gs,ps_data,zkw=0.0)
    @test length(ps_data.zoom_list) == 2
    @test ps_data.zoom_pos == 1

    # OOB: this should cause a warning:
    @test (@test_logs (:warn,r"^unable to zoom") zoomout!(gs,ps_data,zkw=5.0+3.5im)) == -1

    zoomout!(gs,ps_data,zkw=0.5+2.0im)
    @test ps_data.zoom_pos == 1

    @test length(ps_data.zoom_list) == 3
    @test all(map(iscomputed,ps_data.zoom_list))
end
