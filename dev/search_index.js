var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Pseudospectra.jl-1",
    "page": "Home",
    "title": "Pseudospectra.jl",
    "category": "section",
    "text": "Pseudospectra is a Julia package for computing pseudospectra of non-symmetric matrices, and plotting them along with eigenvalues (\"spectral portraits\"). Some related computations and plots are also provided."
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "Whereas the spectrum of a matrix is the set of its eigenvalues, a pseudospectrum is the set of complex numbers \"close\" to the spectrum in some practical sense.More precisely, the ϵ-pseudospectrum of a matrix A, sigma_epsilon(A), is the set of complex numbers lambda such thatlambda is an eigenvalue of some matrix A+E, where the perturbation E is small: E  epsilon\nthe resolvent at lambda has a large norm: (A-λI)^-1  1epsilon,(the definitions are equivalent). Specifically, this package is currently limited to the unweighted 2-norm.Among other things, pseudospectra:elucidate transient behavior hidden to eigen-analysis, and\nindicate the utility of eigenvalues extracted via iterative methods like eigs.This package facilitates computation, display, and investigation of the pseudospectra of matrices and some other representations of linear operators."
},

{
    "location": "index.html#Spectral-portraits-1",
    "page": "Home",
    "title": "Spectral portraits",
    "category": "section",
    "text": "It is customary to display pseudospectra as contour plots of the logarithm of the inverse of the resolvent norm epsilon = 1(A-zI)^-1 for z in a subset of the complex plane. Thus sigma_epsilon(A) is the union of the interiors of such contours. Such plots, sometimes called spectral portraits, are the most prominent product of this package."
},

{
    "location": "index.html#Credit-1",
    "page": "Home",
    "title": "Credit",
    "category": "section",
    "text": "Pseudospectra.jl is largely a translation of the acclaimed MATLAB-based EigTool (homepage here)"
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "The Pseudospectra gateway.\nL.N. Trefethen and M.Embree, Spectra and Pseudospectra; The Behavior of Nonnormal Matrices and Operators, Princeton 2005,"
},

{
    "location": "usage.html#",
    "page": "Usage",
    "title": "Usage",
    "category": "page",
    "text": ""
},

{
    "location": "usage.html#Caveat-1",
    "page": "Usage",
    "title": "Caveat",
    "category": "section",
    "text": "This part of the documentation is a sketch. Please refer to examples and test code."
},

{
    "location": "usage.html#Usage-1",
    "page": "Usage",
    "title": "Usage",
    "category": "section",
    "text": "Typical use of the REPL interface is as follows:using Pseudospectra\nusing Plots\nA = your_matrix_generating_function()\nsetpsplotter()\ngs = setgs()\nps_data = new_matrix(A)\noptions = Dict{Symbol,Any}()\ndriver!(ps_data,options,gs)\n# modify `options` to concentrate on a region of interest\ndriver!(ps_data,options,gs)This should show contour plots of log_10(epsilon) in the vicinity of the spectrum - the standard display of a spectral portrait. Eigenvalues of A are displayed as black points, if they are available.The new_matrix function constructs a stateful object with information about A conducive to pseudospectral analysis; driver! then manages the appropriate computations and plotting."
},

{
    "location": "usage.html#Requirements-1",
    "page": "Usage",
    "title": "Requirements",
    "category": "section",
    "text": "The integrated plotting capabilities require that the Plots and/or PyPlot packages be installed. These are not formal package requirements because much of the Pseudospectra package is useful without plotting."
},

{
    "location": "usage.html#Use-without-the-integrated-plotters-1",
    "page": "Usage",
    "title": "Use without the integrated plotters",
    "category": "section",
    "text": "One can use the psa_compute function to \"simply\" evaluate resolvent norms on a grid, as demonstrated in the plots_minimal.jl script in the examples folder."
},

{
    "location": "lib/public.html#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public.html#Public-interface-1",
    "page": "Public",
    "title": "Public interface",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public.html#Pseudospectra.new_matrix",
    "page": "Public",
    "title": "Pseudospectra.new_matrix",
    "category": "function",
    "text": "new_matrix(A::AbstractMatrix, opts::Dict{Symbol,Any}=()) -> ps_data\n\nprocess a matrix into the auxiliary data structure used by Pseudospectra.\n\nOptions\n\n:direct::Bool: force use of a direct algorithm?\n:keep_sparse::Bool: use sparse matrix code even if A is not large?\n:real_matrix::Bool: treat A as unitarily equivalent to a real matrix?\n:verbosity::Int: obvious\n:eigA: eigenvalues of A, if already known\n:proj_lev: projection level (see psa_compute)\n:npts: edge length of grid for computing and plotting pseudospectra\n:arpack_opts::ArpackOptions: (see type description)\n:levels::Vector{Real}: contour levels\n:ax::Vector{Real}(4): bounding box for computation [xmin,xmax,ymin,ymax]\n:scale_equal::Bool: force isotropic axes for spectral portraits?\n\n\n\n\n\nnew_matrix(A, opts::Dict{Symbol,Any}=()) -> ps_data\n\nprocess a linear operator object into the auxiliary data structure used by Pseudospectra.\n\nThere must be methods with A for eltype, size, and mul!. It is up to the user to make sure that mul! is consistent with any options passed to the iterative solver (see documentation for xeigs).\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Importing-a-matrix-or-operator-1",
    "page": "Public",
    "title": "Importing a matrix or operator",
    "category": "section",
    "text": "new_matrix"
},

{
    "location": "lib/public.html#Pseudospectra.setpsplotter",
    "page": "Public",
    "title": "Pseudospectra.setpsplotter",
    "category": "function",
    "text": "setpsplotter(plotter::Symbol=:default)\n\nSelect a plotting package for use with Pseudospectra.\n\nCurrently :Plots and :PyPlot are implemented. Defaults to :Plots unless PyPlot is already imported without Plots.\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra.setgs",
    "page": "Public",
    "title": "Pseudospectra.setgs",
    "category": "function",
    "text": "setgs(; headless=false, savefigs=true) => gs\n\nConstruct a GUIState for subsequent use by Pseudospectra functions.\n\nAssumes plotting package has been chosen via setpsplotter().\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Setting-up-graphics-subsystem-1",
    "page": "Public",
    "title": "Setting up graphics subsystem",
    "category": "section",
    "text": "setpsplotter\n\nsetgs"
},

{
    "location": "lib/public.html#Pseudospectra.driver!",
    "page": "Public",
    "title": "Pseudospectra.driver!",
    "category": "function",
    "text": "driver!(ps_data, opts, gs; revise_method=false)\n\nCompute pseudospectra and plot a spectral portrait.\n\nIf using an iterative method to get some eigenvalues and a projection, invokes that first.\n\nArguments\n\nps_data::PSAStruct: ingested matrix, as processed by new_matrix\ngs::GUIState: object handling graphical output\nopts::Dict{Symbol,Any}:\n:ax, axis limits (overrides value stored in ps_data).\nother options passed to redrawcontour, arnoldiplotter!\n\nWhen revising a spectral portrait (revise_method==true), the following entries in opts also apply:\n\n:arpack_opts::ArpackOptions,\n:direct::Bool.\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Driver-function-1",
    "page": "Public",
    "title": "Driver function",
    "category": "section",
    "text": "After a matrix has been ingested into a PSAStruct and the graphics subsystem has been established, the following function will compute pseudospectra and plot a spectral portrait:driver!"
},

{
    "location": "lib/public.html#Pseudospectra.psa_compute",
    "page": "Public",
    "title": "Pseudospectra.psa_compute",
    "category": "function",
    "text": "psa_compute(T,npts,ax,eigA,opts,S=I) -> (Z,x,y,levels,info,Tproj,eigAproj,algo)\n\nCompute pseudospectra of a (decomposed) matrix.\n\nUses a modified version of the code in [Trefethen1999]. If the matrix T is upper triangular (e.g. from a Schur decomposition) the solver is much more efficient than otherwise.\n\nArguments\n\nT:      input matrix, usu. from schur()\nnpts:   grid will have npts × npts nodes\nax:     axis on which to plot [min_real, max_real, min_imag, max_imag]\neigA:   eigenvalues of the matrix, usu. also produced by schur(). Pass an empty vector if unknown.\nS:    2nd matrix, if this is a generalized problem arising from an original rectangular matrix.\nopts: a Dict{Symbol,Any} holding options. Keys used here are as follows:\n\nKey Type Default Description\n:levels Vector{Real} auto log10(ϵ) for the desired ϵ levels\n:recompute_levels Bool true automatically recompute ϵ levels?\n:real_matrix Bool eltype(A)<:Real is the original matrix real? (Portrait is symmetric if so.) This is needed because T could be complex even if A was real.\n:proj_lev Real ∞ The proportion by which to extend the axes in all directions before projection. If negative, exclude subspace of eigenvalues smaller than inverse fraction. ∞ means no projection.\n:scale_equal Bool false force the grid to be isotropic?\n\nNotes:\n\nProjection is only done for square, dense matrices.  Projection for sparse matrices may be handled (outside this function) by a Krylov method which reduces the matrix to a projected Hessenberg form before invoking psa_compute.\nThis function does not compute generalized pseudospectra per se. They may be handled by pre- and post-processing.\n\nOutputs:\n\nZ:         the singular values over the grid\nx:         the x coordinates of the grid lines\ny:         the y coordinates of the grid lines\nlevels:   the levels used for the contour plot (if automatically calculated)\nTproj:     the projected matrix (an alias to T if no projection was done)\neigAproj:  eigenvalues projected onto\nalgo: a Symbol indicating which algorithm was used\ninfo:      flag indicating where automatic level creation fails:\n\ninfo Meaning\n0 No error\n-1 No levels in range specified (either manually, or if matrix is too normal to show levels)\n-2 Matrix is so non-normal that only zero singular values were found\n-3 Computation cancelled\n\n[Trefethen1999]: L.N.Trefethen, \"Computation of pseudospectra,\" Acta Numerica 8, 247-295 (1999).\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra-computation-1",
    "page": "Public",
    "title": "Pseudospectra computation",
    "category": "section",
    "text": "psa_compute"
},

{
    "location": "lib/public.html#Pseudospectra.modeplot",
    "page": "Public",
    "title": "Pseudospectra.modeplot",
    "category": "function",
    "text": "modeplot(ps_data,gs,pkey [,z])\n\nExtract and plot an eigenmode or pseudo-eigenmode for the matrix encapsulated in the Pseudospectra object ps_data. Use the value z if provided or prompt for one. If pkey is 1, find the pseudoeigenmode for z; otherwise find the eigenmode for the eigenvalue closest to z.\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Eigen/Pseudo-mode-computation-and-plotting-1",
    "page": "Public",
    "title": "Eigen/Pseudo-mode computation and plotting",
    "category": "section",
    "text": "modeplot"
},

{
    "location": "lib/public.html#Pseudospectra.psa_radius",
    "page": "Public",
    "title": "Pseudospectra.psa_radius",
    "category": "function",
    "text": "psa_radius(A,ϵ [,d]) -> r,z\n\nCompute ϵ-pseudospectral radius for a dense matrix.\n\nQuadratically convergent two-way method to compute the ϵ-pseudospectral radius r of a dense matrix A. Also returns a vector z of points where the pseudospectrum intersects the circle of radius r. Uses the \"criss-cross\" algorithm of Overton and Mengi.\n\nThe ϵ-pseudospectral radius is\n\n   maximum(abs(z)) for z s.t. minimum(σ(A-zI)) == ϵ\n\nOptional arg:\n\nd: eigenvalues of A, if known in advance\n\nKeyword args:\n\n`verbosity: 0 for quiet, 1 for noise\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra.psa_abscissa",
    "page": "Public",
    "title": "Pseudospectra.psa_abscissa",
    "category": "function",
    "text": "psa_abscissa(A,ϵ [,d]) -> α,z\n\nCompute ϵ-pseudospectral abscissa for a dense matrix.\n\nQuadratically convergent two-way method to compute the ϵ-pseudospectral abscissa α of a dense matrix A. Also returns a vector z of points where the pseudospectrum reaches the abscissa. Uses the criss-cross algorithm of Burke et al.\n\nThe ϵ-pseudospectral abscissa is\n\n   maximum(real(z)) for z s.t. minimum(σ(A-zI)) == ϵ\n\nOptional arg:\n\nd: eigenvalues of A, if known in advance\n\nKeyword args:\n\n`verbosity: 0 for quiet, 1 for noise\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra.numerical_range",
    "page": "Public",
    "title": "Pseudospectra.numerical_range",
    "category": "function",
    "text": "numerical_range(A, nstep=20) -> Vector{Complex}\n\nCompute points along the numerical range of a matrix.\n\nNote: this solves an eigensystem for each point, so may be expensive.\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra.numerical_abscissa",
    "page": "Public",
    "title": "Pseudospectra.numerical_abscissa",
    "category": "function",
    "text": "numerical_abscissa(A)\n\nCompute the numerical abscissa of a matrix A, ω(A).\n\nUses eigvals(). ω(A) provides bounds and limiting behavior for norm(expm(t*A)).\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Other-computations-1",
    "page": "Public",
    "title": "Other computations",
    "category": "section",
    "text": "psa_radius\n\npsa_abscissa\n\nnumerical_range\n\nnumerical_abscissa"
},

{
    "location": "lib/public.html#Pseudospectra.surfplot",
    "page": "Public",
    "title": "Pseudospectra.surfplot",
    "category": "function",
    "text": "make a surface plot of the current spectral portrait\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra.mtxexpsplot",
    "page": "Public",
    "title": "Pseudospectra.mtxexpsplot",
    "category": "function",
    "text": "mtxexpsplot(gs::GUIState,ps_data,dt=0.1,nmax=50; gradual=false)\n\nplot the evolution of ∥e^(tA)∥.\n\nThis is useful for analyzing linear initial value problems ∂x/∂t = Ax.\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Pseudospectra.mtxpowersplot",
    "page": "Public",
    "title": "Pseudospectra.mtxpowersplot",
    "category": "function",
    "text": "mtxpowersplot(gs::GUIState,ps_data,nmax=50;gradual=false)\n\nplot norms of powers of a matrix ∥A^k∥\n\nThis is useful for analyzing iterative linear algebra methods.\n\n\n\n\n\n"
},

{
    "location": "lib/public.html#Other-plots-1",
    "page": "Public",
    "title": "Other plots",
    "category": "section",
    "text": "Pseudospectra.surfplot\n\nmtxexpsplot\n\nmtxpowersplot"
},

{
    "location": "lib/internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals.html#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": ""
},

{
    "location": "lib/internals.html#Pseudospectra.Portrait",
    "page": "Internals",
    "title": "Pseudospectra.Portrait",
    "category": "type",
    "text": "Portrait\n\nstructure representing a spectral portrait; includes a mesh for z=x+iy, values of the resolvent norm on the mesh, suitable contour levels, and some other metadata.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#Pseudospectra.PSAStruct",
    "page": "Internals",
    "title": "Pseudospectra.PSAStruct",
    "category": "type",
    "text": "Wrapper structure for Pseudospectra session data\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#Data-structures-1",
    "page": "Internals",
    "title": "Data structures",
    "category": "section",
    "text": "Portrait\n\nPSAStruct"
},

{
    "location": "lib/internals.html#Pseudospectra.xeigs",
    "page": "Internals",
    "title": "Pseudospectra.xeigs",
    "category": "function",
    "text": "xeigs(A, B, channel=nothing; nev=6, ncv=max(20,2*nev+1), which=:LM,\n      tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,)),\n      wantH=true, options=nothing)\n\nCompute (generalized) eigenpairs of A x=λ B x and projections.\n\nModified version of eigs() (from the Arpack package, q.v.) which optionally\n\nprovides the projection matrices,\nprovides intermediate states (typically for plotting).\n\nFor intermediates, caller must provide a Channel argument; in this case xeigs is implemented as a producer which fills channel. When finished, it puts (:finale, d,[v,],nconv,niter,nmult,resid[,H,V]) to channel.\n\nThe type of A must provide methods for mul!, issymmetric, eltype and size. Some variants are not implemented for the case where A is not a factorable matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#Pseudospectra.ArpackOptions",
    "page": "Internals",
    "title": "Pseudospectra.ArpackOptions",
    "category": "type",
    "text": "ArpackOptions{T}(; nev=6, ncv=0, which=:LM,\n                                       tol=zero(T),maxiter=300,\n                                       v0=Vector{T}(0), have_v0=false,\n                                       sigma=nothing)\n\nconstructs an object to manage the Arnoldi scheme; see the documentation for eigs in the Arpack package for the meaning of fields.\n\n\n\n\n\n"
},

{
    "location": "lib/internals.html#Arnoldi-iteration-1",
    "page": "Internals",
    "title": "Arnoldi iteration",
    "category": "section",
    "text": "Pseudospectra.xeigs\n\nArpackOptions"
},

{
    "location": "lib/internals.html#Algorithm-selection-1",
    "page": "Internals",
    "title": "Algorithm selection",
    "category": "section",
    "text": "IRAM means implicitly restarted Arnoldi method, using the ARPACK implementation. \"Large\" and \"small\" are determined by constants which might be made configurable someday."
},

{
    "location": "lib/internals.html#Dense-square-matrix-1",
    "page": "Internals",
    "title": "Dense square matrix",
    "category": "section",
    "text": "If large and not opts[:direct], project via IRAM and go to rectangular branch.\nIf small, use SVD.\nOtherwise, use inverse Lanczos.\nIf a Schur decomposition is available, use that to simplify the solvers."
},

{
    "location": "lib/internals.html#Sparse-square-matrix-1",
    "page": "Internals",
    "title": "Sparse square matrix",
    "category": "section",
    "text": "If small, and not opts[:keep_sparse], convert to dense then use above logic.\nIf opts[:direct], use inverse Lanczos. This involves a sparse factorization for every grid point, so is slow.\nOtherwise project via IRAM and go to rectangular branch."
},

{
    "location": "lib/internals.html#Dense-rectangular-matrix-1",
    "page": "Internals",
    "title": "Dense rectangular matrix",
    "category": "section",
    "text": "If not Hessenberg, factor first.\nIf very tall, factor by QR, otherwise QZ (generalized Schur).\nUse QR for resolvents, then inverse Lanczos.\nUse a simplified QR for large Hessenberg."
},

{
    "location": "lib/demos.html#",
    "page": "Example matrix generators",
    "title": "Example matrix generators",
    "category": "page",
    "text": ""
},

{
    "location": "lib/demos.html#Pseudospectra.landau_fox_li",
    "page": "Example matrix generators",
    "title": "Pseudospectra.landau_fox_li",
    "category": "function",
    "text": "landau_fox_li(N,F)\n\nconstruct a Landau matrix of rank N for Fresnel number F\n\nThe discretization of a Fox-Li operator in the theory of laser cavities by [Landau1977] is one of the earliest applications of pseudospectra.\n\n[Landau1977]: H. J. Landau, \"The notion of approximate eigenvalues applied to an integral equation of laser theory\", Quart. Appl. Math. 35 (1977), pp. 165-172.\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.grcar",
    "page": "Example matrix generators",
    "title": "Pseudospectra.grcar",
    "category": "function",
    "text": "grcar(N)\n\nconstruct a Grcar non-normal matrix of rank N\n\nThis is a popular example in the field of matrix iterations of a matrix whose spectrum is in the right half-plane but whose numerical range is not.  It\'s also a popular example in the study of nonsymmetric Toeplitz matrices.  The matrix was first described in [Grcar1989] and its pseudospectra were first plotted in [Trefethen1991].\n\n[Grcar1989]: J. F. Grcar, \"Operator coefficient methods for linear equations\", tech. report SAND89-8691, Sandia National Labs, 1989\n\n[Trefethen1991]: L. N. Trefethen, \"Psuedospectra of matrices\", in \"Numerical Analysis 1991\" (Dundee 1991), Longman Sci. Tech., Harlow, 1992, 234-266.\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.demmel",
    "page": "Example matrix generators",
    "title": "Pseudospectra.demmel",
    "category": "function",
    "text": "demmel(N,M)\n\nconstruct Demmel\'s matrix of rank N with edge value M\n\nDemmel\'s matrix was perhaps the first matrix whose pseudospectra (with N=3) appeared in [Demmel1987] Demmel devised the matrix in order to disprove a conjecture of Van Loan\'s [VanLoan1985].\n\n[Demmel1987]: J. W. Demmel, \"A counterexample for two conjectures about stability\", IEEE Trans. Auto. Control, AC-32, pp. 340-343, 1987.\n\n[VanLoan1985]: C. Van Loan, \"How near is a stable matrix to an unstable matrix?\", in R. Brualdi, et al., eds, Contemporary Mathematics 47, Amer. Math. Soc., 1985\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.supg",
    "page": "Example matrix generators",
    "title": "Pseudospectra.supg",
    "category": "function",
    "text": "supg(N)\n\nconstruct a SUPG matrix for an NxN mesh\n\nThe matrix arises from SUPG discretisation of an advection-diffusion operator.\n\nThis demo is based on a function by Mark Embree, which was essentially distilled from the IFISS software of Elman and Silvester.  See [Fischer1999] and [DJS].\n\nFor this example, whilst the eigenvalues computed by eigs are inaccurate due to the sensitivity of the problem, the approximate pseudospectra are very good approximations to the true ones.\n\n[Fischer1999]: B. Fischer, A. Ramage, D. Silvester and A. Wathen, \"Towards Parameter-Free Streamline Upwinding for Advection-Diffusion Problems\", Comput. Methods Appl. Mech. Eng., 179, 1999, 185-202.\n\n[DJS]: http://www.ma.umist.ac.uk/djs/software.html\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.olmstead",
    "page": "Example matrix generators",
    "title": "Pseudospectra.olmstead",
    "category": "function",
    "text": "olmstead(N)\n\nload Olmstead matrix from file (N ∈ {500,1000})\n\nnote: Note\nThis function requires a Matrix Market file from [http://math.nist.gov/MatrixMarket/data/NEP/olmstead/olmstead.html]\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.convdiff_fd",
    "page": "Example matrix generators",
    "title": "Pseudospectra.convdiff_fd",
    "category": "function",
    "text": "convdiff_fd(N)\n\nconstruct matrix representation of 2-D NxN convection-diffusion operator\n\nconvdiff_fd(N) returns matrix C of dimension N^2 representing the operator          -ϵ ∇^2 u  + ∂u/∂y\n\non the unit square with Dirichlet boundary conditions. The discretisation is done using finite differences (see, e.g., [Hemmingsson1998]), resulting in large, sparse matrices. The code is based on a routine by Daniel Loghin.\n\n[Hemmingsson1998]: L. Hemmingsson, \"A Semi-circulant Preconditioner for the Convection-Diffusion Equation\", Numer. Math 81, 211-248 (1998).\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.advdiff",
    "page": "Example matrix generators",
    "title": "Pseudospectra.advdiff",
    "category": "function",
    "text": "advdiff(N,η=0.015)\n\nconstruct the matrix representing an advection-diffusion operator via Chebyshev collocation.\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.orrsommerfeld",
    "page": "Example matrix generators",
    "title": "Pseudospectra.orrsommerfeld",
    "category": "function",
    "text": "orrsommerfeld(N,R=5722,α=1.0)\n\nconstruct a matrix representing the Orr-Sommerfeld operator in the rank N+1 Chebyshev value-space basis. Note: this imposes a norm of dubious physical meaning, since a clever trick is used to suppress spurious eigenvalues.\n\n[HS1994] W.Z.Huang and D.M.Sloan, \"The pseudospectral method for solving differential eigenvalue equations,\" J. Comput. Phys. 111, 399-409 (1994).\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.rdbrusselator",
    "page": "Example matrix generators",
    "title": "Pseudospectra.rdbrusselator",
    "category": "function",
    "text": "rdbrusselator(N)\n\nload reaction-diffusion Brusselator matrix from file (N ∈ {800,3200})\n\nnote: Note\nThis function requires a Matrix Market file from [http://math.nist.gov/MatrixMarket/data/NEP/brussel/brussel.html]\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.trefethen_tutorial",
    "page": "Example matrix generators",
    "title": "Pseudospectra.trefethen_tutorial",
    "category": "function",
    "text": "trefethen_tutorial(N) -> B,A\n\ncompute the matrices used for tutorial purposes in [Trefethen1999]. A is the Chebyshev discretization of a Schrödinger operator with complex potential. B includes weight factors so that a basic L2 norm is appropriate.\n\n[Trefethen1999]: L.N.Trefethen, \"Computation of Pseudospectra\", Acta Numerica 1999.\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Pseudospectra.kahan",
    "page": "Example matrix generators",
    "title": "Pseudospectra.kahan",
    "category": "function",
    "text": "kahan(N)\n\nconstruct Kahan\'s matrix of rank N\n\nKahan\'s matrix [Kahan1966] was devised to illustrate that QR factorisation with column pivoting is not a fail-safe method for determining the rank of the matrix. Rank determination is related to the \"distance to singularity\" of a matrix [Higham1988], or in other words, how large a perturbation is needed to make the matrix singular. The pseudospectra are shown in [Trefethen1991].\n\n[Higham1988]: N. J. Higham, \"Matrix nearness problems and applications\", NA Rep. 161, Dept. of Maths., U. of Manchester, June 1988.\n\n[Kahan1966]: W. Kahan, \"Numerical linear algebra\", Canad. Math. Bull. 9, pp. 757-801, 1966.\n\n[Trefethen1991]: L. N. Trefethen, \"Psuedospectra of matrices\", in \"Numerical Analysis 1991\" (Dundee 1991), Longman Sci. Tech., Harlow, 1992, 234-266.\n\n\n\n\n\n"
},

{
    "location": "lib/demos.html#Example-matrix-generators-1",
    "page": "Example matrix generators",
    "title": "Example matrix generators",
    "category": "section",
    "text": "These are mainly translations of demonstrations from EigTool.note: Note\nThese functions are not exported.Pseudospectra.landau_fox_li\n\nPseudospectra.grcar\n\nPseudospectra.demmel\n\nPseudospectra.supg\n\nPseudospectra.olmstead\n\nPseudospectra.convdiff_fd\n\nPseudospectra.advdiff\n\nPseudospectra.orrsommerfeld\n\nPseudospectra.rdbrusselator\n\nPseudospectra.trefethen_tutorial\n\nPseudospectra.kahan"
},

]}
