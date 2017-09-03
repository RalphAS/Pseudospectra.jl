#=
This is a template for a development/debugging script.
=#

# Wrap things up in a module so we can reload stuff as needed.
module M

using Pseudospectra, Base.Test

# To run like a test:
# displayplots = true # If you want to see results immediately.
include("../test/testplotter.jl")

# Alternatively, to run like a user:
#=
using PyPlot
setpsplotter(:PyPlot)
gs = setgs()
=#

A = # something new and challenging

opts = Dict{Symbol,Any}(:verbosity => 2,
                        # whatever
                        )
ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

end # module
