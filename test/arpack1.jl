module M
# 400x400 real sparse iterative

using Pseudospectra
# using Base.Test

include("testplotter.jl")

A = Pseudospectra.convdiff_fd(10)
opts = Dict{Symbol,Any}(:ARPACK_plots => false,:npts=>40)
opts[:arpack_opts] = ArpackOptions{eltype(A)}(which=:SR,nev=20,ncv=60)

ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

end
