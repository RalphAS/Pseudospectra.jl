#=
Set up the infrastructure
=#
Pkg.add("Plots")
using Plots
pyplot()
p=scatter(rand(10),rand(10))
png(p,"/tmp/fubar.png")
