#=
This is a template for a debugging script.
=#

# Remember to set displayplots = true if you want to see results.

# Wrap things up in a module so we can reload stuff as needed.
module M

using Pseudospectra
include("../test/testplotter.jl")

A = # something new and challenging


opts = Dict{Symbol,Any}(:verbosity => 2,
                        # whatever
                        )
ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

end # module
