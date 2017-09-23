"""
object to hold state for the GUI used by Pseudospectra
"""
@compat abstract type GUIState end

type LevelDesc
    isunif::Bool
    full_levels
    first
    step
    last
    function LevelDesc(levs)
        if isa(levs,Range)
            levels = collect(levs)
        elseif isa(levs,Vector) && (eltype(levs) <: Real)
            levels = levs
        else
            throw(ArgumentError("level specification must be a range or real vector"))
        end
        isunif = ((length(levels) == 2)
                  || ((length(levels) > 2)
                      && (maximum(abs.(diff(diff(levels)))) < 10*eps())))
        if isunif
            first = levels[1]
            step = levels[2]-levels[1]
            last = levels[end]
            full_levels = nothing
        else
            full_levels = copy(levels)
            first = nothing; step = nothing; last = nothing
        end
        new(isunif,full_levels,first,step,last)
    end
end

"""
    Portrait

structure representing a spectral portrait; includes a mesh for `z=x+iy`,
values of the resolvent norm on the mesh, suitable contour levels,
and some other metadata.
"""
type Portrait
    x::Vector
    y::Vector
    Z::Matrix
    npts::Int
    ax::Vector
    levels::LevelDesc
    autolev::Bool
    proj_lev
    dims
    computed::Bool
    scale_equal::Bool
end

import Base: show

function show(io::IO,z::Portrait)
    k = isempty(z.Z) ? "missing" : "present"
    print(io,"npts: $(z.npts), ax: $(z.ax), autolev: $(z.autolev), "
          * "computed: $(z.computed), data $k, proj $(z.proj_lev), "
          * "scale_eq: $(z.scale_equal)\n",
          "levels: $(z.levels)")
end

"""
    ArpackOptions{T}(; nev=6, ncv=0, which=:LM,
                                           tol=zero(T),maxiter=300,
                                           v0=Vector{T}(0), have_v0=false,
                                           sigma=nothing)

constructs an object to manage the Arnoldi scheme; see [`eigs`](@ref) for the
meaning of fields.
"""
type ArpackOptions{T}
    # Control of iterative computations:
    nev::Int # nr. of eigenvalues for eigs() to search for
    ncv::Int # max. subspace size for eigs()
    which::Symbol # as for eigs()
    tol::Real # tolerance for eigs()
    maxiter::Int # bound for eigs()
    v0::Vector # initial Ritz vector for eigs()
    sigma # shift parameter for eigs()

    function (::Type{ArpackOptions{T}}){T}(; nev=6, ncv=0, which=:LM,
                                           tol=zero(T),maxiter=300,
                                           v0=Vector{T}(0),
                                           sigma=nothing)
        if ncv==0
            ncv = max(20,2*nev+1)
        end
        new{T}(nev,ncv,which,tol,maxiter,v0,sigma)
    end
end

function Base.:(==){T,S}(l::ArpackOptions{T},r::ArpackOptions{S})
    l === r && return true
    for f in fieldnames(ArpackOptions)
        (getfield(l,f) == getfield(r,f)) || return false
    end
    true
end
function Base.hash{T}(l::ArpackOptions{T},h::UInt)
    h1 = hash(:ArpackOptions,h)
    for f in fieldnames(ArpackOptions)
        h1 = hash(getfield(l,f),h1)
    end
    h1
end

"""
Wrapper structure for Pseudospectra session data
"""
type PSAStruct
    matrix
    unitary_mtx
    input_matrix
    input_unitary_mtx
    ps_dict::Dict{Symbol,Any}
    zoom_list::Vector{Portrait}
    zoom_pos::Int
    # aside: why is there a :proj_lev entry in ps_dict?
    function PSAStruct(m1,u1,m1i,u1i,dict)
        new(m1,u1,m1i,u1i,dict,Vector{Portrait}(0),0)
    end
end
#=
 A (probably incomplete) list of keys in ps_dict
 :Aisreal
 :isHessenberg
 :schur_mtx
 :schur_unitary_mtx
 :direct
 :ews
 :orig_ews
 :ew_estimates
 :arpack_opts
 :init_opts
 :init_direct
 :verbosity
 :fov - samples of the numerical range (field of values)
 :projection_on
 :proj_lev
 :proj_valid
 :proj_ews
 :proj_matrix - the Hessenberg matrix of coeffts from IRAM (ARPACK)
 :proj_unitary_mtx - the Krylov basis from IRAM (ARPACK)
 :matrix2 - second matrix (if QZ is used for rectangular A)

 :reentering - flag for logic (not in EigTool)
 :proj_axes
 :comp_proj_lev
 :mode_markers

 Used only in GUI code:
 :zpsradius
 :zpsabscissa
=#
