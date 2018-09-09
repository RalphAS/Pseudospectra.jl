# pseudospectral radius and abscissa
using Pseudospectra, Test

@testset "PS radius & abscissa" begin
    A = Pseudospectra.grcar(100)

    opts = Dict{Symbol,Any}()

    # also check that this gets the eigenvalues
    ps_data = new_matrix(A,opts)
    topeig = maximum(imag(ps_data.ps_dict[:ews]))
    @test isapprox(topeig, 2.2617668, atol=1e-6)

    # reference values from Guglielmi & Overton (2012)
    r,zr = psa_radius(A,1e-4,ps_data.ps_dict[:ews])
    @test isapprox(r, 2.85216, atol=1e-5)
    α,za = psa_abscissa(A,1e-4,ps_data.ps_dict[:ews])
    @test isapprox(α, 2.41276, atol=1e-5)

    # Note: if we ever interface to a Hamiltonian eigensolver, we could use
    # the ϵ=1e-8 case in S&P.
end
