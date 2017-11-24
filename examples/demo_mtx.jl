#=
This file is part of Pseudospectra.jl.

Julia translation
copyright 2017 Ralph Smith

Portions derived from EigTool
Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
of the University of Oxford, and the EigTool Developers. All rights reserved.

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#

"""
    grcar(N)

construct a Grcar non-normal matrix of rank N

This is a popular example in the field of matrix iterations
of a matrix whose spectrum is in the right half-plane but
whose numerical range is not.  It's also a popular example
in the study of nonsymmetric Toeplitz matrices.  The matrix
was first described in [^Grcar1989] and its pseudospectra were first
plotted in [^Trefethen1991].

[^Grcar1989]: J. F. Grcar, "Operator coefficient methods for linear equations", tech. report SAND89-8691, Sandia National Labs, 1989
[^Trefethen1991]: L. N. Trefethen, "Psuedospectra of matrices", in "Numerical Analysis 1991" (Dundee 1991), Longman Sci. Tech., Harlow, 1992, 234-266.
"""
function grcar(N::Integer)
  G = (diagm(ones(N),0) - diagm(ones(N-1),-1) +
      diagm(ones(N-1),1) + diagm(ones(N-2),2) +
      diagm(ones(N-3),3))
end

"""
    cheb(N)

construct a matrix representing differentiation of rank N+1 Chebyshev
value-space basis, based on [Trefethen2000]

[^Trefethen2000] L.N.Trefethen, *Spectral Methods in MATLAB*, SIAM, Philadelphia, 2000.
"""
function cheb(N::Integer)
    if N==0
        return 0.0,1.0
    end
    x = cos.((π/N)*collect(0:N))
    c = vcat(2,ones(N-1),2) .* map(x -> (-1)^x,0:N)
    dX = x .- x'
    D = (c*(1 ./ c)') ./ (dX + eye(N+1))
    D = D - diagm(squeeze(sum(D,2),2))
    return D,x
end

"""
    advdiff(N,η=0.015)

construct the matrix representing an advection-diffusion operator via
Chebyshev collocation.
"""
function advdiff(N,η=0.015)
    D,x = cheb(N)
    # rescale to interval 0:1
    D = 2*D
    x = (x+1)*(1/2)
    L = η*D^2 + D
    # impose boundary conditions
    L = L[2:N,2:N]
    # Gauss-Chebyshev quadratrure weights
    w = sqrt.(π*sqrt.(x-x.^2)*(1/(2N)))
    # convert to L2 norm on function space
    w = w[2:N]
    L = diagm(w) * L * diagm(1.0 ./ w)
end

"""
    orrsommerfeld(N,R=5722,α=1.0)

construct a matrix representing the Orr-Sommerfeld operator in the
rank N+1 Chebyshev value-space basis.
Note: this imposes a norm of dubious physical meaning, since a clever
trick is used to suppress spurious eigenvalues.

[^HS1994] W.Z.Huang and D.M.Sloan, "The pseudospectral method for solving differential eigenvalue equations," J. Comput. Phys. 111, 399-409 (1994).
"""
function orrsommerfeld(N,R=5722,α=1.0)
    D,x = cheb(N)
    D2 = D^2
    D2 = D2[2:N,2:N]
    S = diagm([0; 1 ./(1-x[2:N].^2); 0])
    D4 = (diagm(1-x.^2)*D^4 - 8*diagm(x)*D^3 - 12*D^2)*S
    D4 = D4[2:N,2:N]
    a2 = α^2
    I = eye(N-1)
    A = (D4-2*a2*D2+a2^2*I)*(1/R) - 2α*im*I - α*im*diagm(1-x[2:N].^2)*(D2-a2*I)
    B = D2-a2*I
    C = B\A
end

"""
    trefethen_tutorial(N) -> B,A

compute the matrices used for tutorial purposes
in [^Trefethen1999]. `A` is the Chebyshev discretization
of a Schrödinger operator with complex potential.
`B` includes weight factors so that a basic L2 norm is appropriate.

[^Trefethen1999]: L.N.Trefethen, "Computation of Pseudospectra", Acta Numerica 1999.
"""
function trefethen_tutorial(N)
    D = zeros(N+2,N+2)
    iran = 0:(N+1)
    ci = vcat(2,ones(N),2)
    x = cos(π*iran*(1/(N+1)))
    for j=0:N+1
        cj = ((j==0) || (j==N+1)) ? 2 : 1
        denom = cj * (x[iran+1]-x[j+1])
        denom[j+1] = 1
        D[iran+1,j+1] = ci .* (-1.0).^(iran+j) ./ denom
        if (j > 0) && (j<N+1)
            D[j+1,j+1] = -(1/2)*x[j+1] / (1-x[j+1]^2)
        end
    end
    D[1,1] = (2*(N+1)^2+1)/6
    D[N+2,N+2] = -(2*(N+1)^2+1)/6
    L = 10.0
    D = D / L
    A = D^2
    # absorb boundary conditions
    x = x[2:N+1]
    x = L*x
    A = A[2:N+1,2:N+1]
    A = A + (3+3im) * diagm(x.^2) - (1/16)*diagm(x.^4)
    # absorb Chebyshev weights
    w = sqrt.(π*sqrt(L^2 - x.^2) * (1 / (2*(N+1))))
    B = zeros(eltype(A),N,N)
    for j=1:N
        B[:,j] = w .* A[:,j] * (1/w[j])
    end
    B,A
end


"""
    convdiff_fd(N)

construct matrix representation of 2-D NxN convection-diffusion operator

convdiff_fd(N) returns matrix C of dimension N^2 representing the
operator
         -ϵ ∇^2 u  + ∂u/∂y

on the unit square with Dirichlet boundary conditions. The
discretisation is done using finite differences (see, e.g., [^Hemmingsson1998]),
resulting in large, sparse matrices. The code is based on a
routine by Daniel Loghin.

[^Hemmingsson1998]: L. Hemmingsson, "A Semi-circulant Preconditioner for the Convection-Diffusion Equation", Numer. Math 81, 211-248 (1998).

"""
function convdiff_fd(N)

    # Diffusion parameter
    e = 1e-2

    # Wind blowing in y-direction.
    w1 = 0
    w2 = 1

    # Set up the grid size (more points in y direction due to wind)
    n1 = N
    n2 = 4*N
    h1 = 1/(n1+1)
    h2 = 1/(n2+1)

    # Generate 1D finite difference matrices
    alpha1 = w1/h1
    beta1 = 2*e/h1^2
    alpha2 = w2/h2
    beta2 = 2*e/h2^2

    A1 = spdiagm(((-alpha1-beta1)*ones(n1-1),
                  2*beta1*ones(n1),
                  (alpha1-beta1)*ones(n1-1)), -1:1, n1, n1)
    A2 = spdiagm(((-alpha2-beta2)*ones(n2-1),
                  2*beta2*ones(n2),
                  (alpha2-beta2)*ones(n2-1)), -1:1, n2, n2)

    # Now compute the 2D finite difference matrix
    C = kron(speye(n2),A1) + kron(A2,speye(n1))
end

function gallery5(T=Float64)
    @assert T <: Real
    A = hcat([-9,70,-575,3891,1024],
             [11,-69,575,-3891,-1024],
             [-21,141,-1149,7782,2048],
             [63,-421,3451,-23345,-6144],
             [-252,1684,-13801,93365,24572])
    Matrix{T}(A)
end

"""
    kahan(N)

construct Kahan's matrix of rank N

Kahan's matrix [^Kahan1966] was devised to illustrate that QR
factorisation with column pivoting is not a fail-safe
method for determining the rank of the matrix. Rank
determination is related to the "distance to singularity"
of a matrix [^Higham1988], or in other words, how large a perturbation
is needed to make the matrix singular. The pseudospectra are
shown in [^Trefethen1991].

[^Higham1988]: N. J. Higham, "Matrix nearness problems and applications", NA Rep. 161, Dept. of Maths., U. of Manchester, June 1988.
[^Kahan1966]: W. Kahan, "Numerical linear algebra", Canad. Math. Bull. 9, pp. 757-801, 1966.
[^Trefethen1991]: L. N. Trefethen, "Psuedospectra of matrices", in "Numerical Analysis 1991" (Dundee 1991), Longman Sci. Tech., Harlow, 1992, 234-266.
"""
function kahan(N)
    s = (1/10)^(1/(N-1))
    c = sqrt(1-s^2)

    col = -s.^(0:N-1)*c
    K = triu(repmat(col,1,N),1)+diagm(s.^(0:N-1))
end

"""
    demmel(N,M)

construct Demmel's matrix of rank N with edge value M

Demmel's matrix was perhaps the first matrix
whose pseudospectra (with N=3) appeared in [^Demmel1987]
Demmel devised the matrix in order to disprove
a conjecture of Van Loan's [^VanLoan1985].

[^Demmel1987]: J. W. Demmel, "A counterexample for two conjectures about stability", IEEE Trans. Auto. Control, AC-32, pp. 340-343, 1987.
[^VanLoan1985]: C. Van Loan, "How near is a stable matrix to an unstable matrix?", in R. Brualdi, et al., eds, Contemporary Mathematics 47, Amer. Math. Soc., 1985
"""
function demmel(N,M=1e4)
    isdefined(Main, :ToeplitzMatrices) ||
        throw(ErrorException("this function requires the ToeplitzMatrices package"))
  B = M^(1/(N-1))
  c = vcat(ones(1), zeros(N-1))
  r = B.^(0:N-1)
  D = -full(Main.ToeplitzMatrices.Toeplitz(c,r))
end

"""
    landau_fox_li(N,F)

construct a Landau matrix of rank N for Fresnel number F

The discretization of a Fox-Li operator in the theory of laser cavities
by [^Landau1977] is one of the earliest applications of pseudospectra.

[^Landau1977]: H. J. Landau, "The notion of approximate eigenvalues applied to an integral equation of laser theory", Quart. Appl. Math. 35 (1977), pp. 165-172.
"""
function landau_fox_li(N,F=0)
    (F == 0) && (F = (N > 200) ? 32 : 12)
    # calculate Gaussian quadrature nodes and weights
    β = 0.5 * (1 - (2.0 * collect(1:N-1)).^(-2)) .^(-1/2)
    T = diagm(β,1) + diagm(β,-1)
    d,V = eig(T)
    idx = sortperm(d)
    nodes = d[idx]
    w = zeros(N)
    w[1:N] = 2*V[1,idx].^2
    B = zeros(N,N)+0im
    for k in 1:N
        B[k,:] = w' * sqrt(F*im) .* exp.(-im*π*F*(nodes[k] - nodes').^2)
    end
    w = sqrt.(w)
    for j in 1:N
        B[:,j] = w .* B[:,j] * (1/w[j])
    end
    B
end

"""
    supg(N)

construct a SUPG matrix for an NxN mesh

The matrix arises from SUPG discretisation of an advection-diffusion operator.

This demo is based on a function by Mark Embree, which
was essentially distilled from the IFISS software of
Elman and Silvester.  See [^Fischer1999] and [^DJS].

For this example, whilst the eigenvalues computed by eigs
are inaccurate due to the sensitivity of the problem, the
approximate pseudospectra are very good approximations to
the true ones.

[^Fischer1999]: B. Fischer, A. Ramage, D. Silvester and A. Wathen, "Towards Parameter-Free Streamline Upwinding for Advection-Diffusion Problems", Comput. Methods Appl. Mech. Eng., 179, 1999, 185-202.
[^DJS]: http://www.ma.umist.ac.uk/djs/software.html
"""
function supg(N)
    h = 1/(N+1)
    theta = 0
    nu = (N < 20) ? 0.01 : 0.0001

    # Compute the "optimal" upwinding parameter
    hbar = h/(max(sin(theta), cos(theta)))
    delta = max(0, 1/2 - nu/hbar)

    w = [-sin(theta),cos(theta)]
    wx2 = w[1]^2
    wy2 = w[2]^2

    if (N > 10)
        A = spzeros(N*N, N*N)
        v = ones(N)
        for j=1:N
            blk =  v*(nu*[-1/3 8/3 -1/3] +
                      h*[-w[1]/3  0 w[1]/3] +
                      delta*h*[(wy2-2*wx2)/3 4/3*(wx2+wy2) (wy2-2*wx2)/3])
            A[N*(j-1)+1:N*j, N*(j-1)+1:N*j] = spdiagm((blk[1:end-1,1],
                                                       blk[:,2],
                                                       blk[1:end-1,3]),
                                                      [-1,0,1],N,N)
            if (j>1)
                blk = v*(nu*[-1/3 -1/3 -1/3] +
                         h*[-(w[1]+w[2])/12 -w[2]/3 (w[1]-w[2])/12] +
                         delta*h*
                         [(-(wx2+wy2)/6-w[1]*w[2]/2) (wx2-2*wy2)/3 (-(wx2+wy2)/6+w[1]*w[2]/2)])
                A[N*(j-1)+1:N*j, N*(j-2)+1:N*(j-1)] = spdiagm((blk[1:end-1,1],
                                                               blk[:,2],
                                                               blk[1:end-1,3]),
                                                              [-1,0,1],N,N)
            end
            if (j<N)
                blk = v*(nu*[-1/3 -1/3 -1/3] +
                         h*[(w[2]-w[1])/12 w[2]/3 (w[1]+w[2])/12] +
                         delta*h*
                         [(-(wx2+wy2)/6+w[1]*w[2]/2) (wx2-2*wy2)/3 (-(wx2+wy2)/6-w[1]*w[2]/2)])
                A[N*(j-1)+1:N*j, N*j+1:N*(j+1)] = spdiagm((blk[1:end-1,1],
                                                           blk[:,2],
                                                           blk[1:end-1,3]),
                                                          [-1,0,1],N,N)
            end
        end

    else
        nn = N*N
        H = zeros(nn,nn)
        S = zeros(nn,nn)
        U = zeros(nn,nn)

        for j=1:N
            for k=1:N
                row = (j-1)*N+k
                H[row, row] = 8/3
                (k > 1) && (H[row, row-1] = -1/3)
                (k < N) && (H[row, row+1] = -1/3)
                (j > 1) && (H[row, row-N] = -1/3)
                (j < N) && (H[row, row+N] = -1/3)
                ((j > 1) & (k > 1)) && (H[row, row-N-1] = -1/3)
                ((j > 1) & (k < N)) && (H[row, row-N+1] = -1/3)
                ((j < N) & (k > 1)) && (H[row, row+N-1] = -1/3)
                ((j < N) & (k < N)) && (H[row, row+N+1] = -1/3)
            end
        end

        for j=1:N
            for k=1:N
                row = (j-1)*N+k
                S[row, row] = 0
                (k > 1) && (S[row, row-1] = -1/3*w[1])
                (k < N) && (S[row, row+1] =  1/3*w[1])
                (j > 1) && (S[row, row-N] = -1/3*w[2])
                (j < N) && (S[row, row+N] =  1/3*w[2])
                ((j > 1) & (k > 1)) && (S[row, row-N-1] = -(w[1]+w[2])/12)
                ((j > 1) & (k < N)) && (S[row, row-N+1] =  (w[1]-w[2])/12)
                ((j < N) & (k > 1)) && (S[row, row+N-1] =  (w[2]-w[1])/12)
                ((j < N) & (k < N)) && (S[row, row+N+1] =  (w[1]+w[2])/12)
            end
        end

        for j=1:N
            for k=1:N
                row = (j-1)*N+k
                U[row, row] = 4/3*(wx2+wy2)
                (k > 1) && (U[row, row-1] = (wy2-2*wx2)/3)
                (k < N) && (U[row, row+1] = (wy2-2*wx2)/3)
                (j > 1) && (U[row, row-N] = (wx2-2*wy2)/3)
                (j < N) && (U[row, row+N] = (wx2-2*wy2)/3)
                ((j > 1) & (k > 1)) && (U[row, row-N-1] = -(wx2+wy2)/6 - w[1]*w[2]/2)
                ((j > 1) & (k < N)) && (U[row, row-N+1] = -(wx2+wy2)/6 + w[1]*w[2]/2)
                ((j < N) & (k > 1)) && (U[row, row+N-1] = -(wx2+wy2)/6 + w[1]*w[2]/2)
                ((j < N) & (k < N)) && (U[row, row+N+1] = -(wx2+wy2)/6 - w[1]*w[2]/2)
            end
        end

        A = nu*H + h*S + delta*h*U
    end
    A
end # supg

"""
    olmstead(N)

load Olmstead matrix from file (`N ∈ {500,1000}`)

!!! note

    This function requires a Matrix Market file from [http://math.nist.gov/MatrixMarket/data/NEP/olmstead/olmstead.html]
"""
function olmstead(N)
    isdefined(Main, :MatrixMarket) ||
        throw(ErrorException("this function requires the MatrixMarket package"))
    if N==500
        A = Main.MatrixMarket.mmread(joinpath(dirname(@__FILE__),"data","olm500.mtx"))
    elseif N==1000
        A = Main.MatrixMarket.mmread(joinpath(dirname(@__FILE__),"data","olm1000.mtx"))
    else
        throw(ArgumentError("no data available for size $N"))
    end
    A
end
"""
    rdbrusselator(N)

load reaction-diffusion Brusselator matrix from file (`N ∈ {800,3200}`)

!!! note

    This function requires a Matrix Market file from [http://math.nist.gov/MatrixMarket/data/NEP/brussel/brussel.html]
"""
function rdbrusselator(N)
    isdefined(Main, :MatrixMarket) ||
        throw(ErrorException("this function requires the MatrixMarket package"))
    if N==800
        A = Main.MatrixMarket.mmread(joinpath(dirname(@__FILE__),"data","rdb800l.mtx"))
    elseif N==3200
        A = Main.MatrixMarket.mmread(joinpath(dirname(@__FILE__),"data","rdb3200l.mtx"))
    else
        throw(ArgumentError("no data available for size $N"))
    end
    A
end
