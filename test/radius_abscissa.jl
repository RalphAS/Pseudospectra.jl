# pseudospectral radius and abscissa

# Note: if we ever interface to a Hamiltonian eigensolver, we could use
# the ϵ=1e-8 case in S&P.

using Pseudospectra, Test
using Pseudospectra: selectfig

@testset "PS radius & abscissa" begin
    A = Pseudospectra.grcar(100)

    opts = Dict{Symbol,Any}()

    have_canvas = psplotter ∈ [:PythonPlot, :PyPlot]
    have_canvas && selectfig(gs, false)
    plotfig = have_canvas ? gs.secondaryfignum : nothing
    # also check that this gets the eigenvalues
    ps_data = new_matrix(A,opts)
    topeig = maximum(imag(ps_data.ps_dict[:ews]))
    @test isapprox(topeig, 2.2617668, atol=1e-6)

    # reference values from Guglielmi & Overton (2012)
    r,zr = psa_radius(A,1e-4,ps_data.ps_dict[:ews]; fig_id=plotfig)
    have_canvas && gs.drawcmd(gs, gs.secondaryph, 1)
    @test isapprox(r, 2.85216, atol=1e-5)
    α,za = psa_abscissa(A,1e-4,ps_data.ps_dict[:ews]; fig_id=plotfig)
    have_canvas && gs.drawcmd(gs, gs.secondaryph, 1)
    @test isapprox(α, 2.41276, atol=1e-5)

    have_canvas && selectfig(gs, true)
end
