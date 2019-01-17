#=
Set up the infrastructure, and trigger a quick failure if appropriate.
=#
if Sys.isapple()
    using Pkg
    Pkg.add("Conda")
    using Conda
    Conda.add("pyqt")
end
Pkg.add("Plots")
using Plots
pyplot()
p=scatter(rand(10),rand(10))
png(p,"/tmp/fubar.png")
