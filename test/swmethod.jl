using Pseudospectra, Test

@testset "Switch-method" begin
    # 400x400 real sparse
    A = Pseudospectra.convdiff_fd(10)



    # first direct, then ARPACK
    opts = Dict{Symbol,Any}(:direct => true,:npts=>24,:ax => [-10,100,-50,50])
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    ps_data_dir = deepcopy(ps_data)

    Pseudospectra.set_method!(ps_data,false)
    opts = Dict{Symbol,Any}(:direct => false, :ARPACK_plots => false,:npts=>40,
                            :ax => [30,100,-50,50])
    opts[:arpack_opts] = ArpackOptions{eltype(A)}(which=:SR,nev=20,ncv=60)

    driver!(ps_data,opts,gs,revise_method=true)
    @test iscomputed(ps_data)

    # first ARPACK, then direct
    opts = Dict{Symbol,Any}(:ARPACK_plots => false,:npts=>40,
                            :ax => [30,100,-50,50])
    opts[:arpack_opts] = ArpackOptions{eltype(A)}(which=:SR,nev=20,ncv=60)
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)
    @test iscomputed(ps_data)
    ps_data_arpack = deepcopy(ps_data)

    Pseudospectra.set_method!(ps_data,true)

    opts = Dict{Symbol,Any}(:direct => true,:npts=>24,:ax => [-10,100,-50,50])
    driver!(ps_data,opts,gs,revise_method=true)
    @test iscomputed(ps_data)

    # first ARPACK, then with new opts
    ps_data = deepcopy(ps_data_arpack)

    opts = Dict{Symbol,Any}(:ARPACK_plots => false,
                            :npts=>24,:ax => [-10,100,-50,50])
    opts[:arpack_opts] = ArpackOptions{eltype(A)}(which=:SR,nev=40,ncv=80)

    driver!(ps_data,opts,gs,revise_method=true)
    @test iscomputed(ps_data)
    # CHECKME: if 41 are computed, are all saved?
    @test length(ps_data.ps_dict[:ews]) >= 40
end
