#=
This file is part of Pseudospectra.jl.
See LICENSE file for details.

Julia translation
copyright (c) 2017 Ralph Smith

Portions derived from EigTool
Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
of the University of Oxford, and the EigTool Developers. All rights reserved.
=#

"""
    recalc_levels(σ, ax)

construct a reasonable set of contour levels for displaying `log10(σ)`
where `σ` is a meshed field on axes `ax`.
"""
function recalc_levels(sigmin,ax)
    err = 0
    smin,smax = extrema(sigmin)
    if smax <= smallσ
        levels = Vector{typeof(smin)}(0)
        err = -2
        return levels,err
    end
    if smin <= smallσ
        smin=minimum(sigmin[sigmin .>= smallσ])
    end
    max_val = log10(smax)
    min_val = log10(smin)
    num_lines = 8
    scale = min(ax[4]-ax[3],ax[2]-ax[1])
    last_lev = log10(0.03*scale)
    if max_val >= last_lev
        num_lines = ceil(Int,num_lines*(last_lev-min_val)/(max_val-min_val))
        max_val = last_lev
    end
    if max_val < min_val
        max_val = log10(smin) + 0.1*log10(smax/smin)
        num_lines = 3
    end
    max_lines = num_lines
    if num_lines > 0
        stepsize = max(round((max_val-min_val)/num_lines*4)/4,0.25)
        (stepsize >= 0.75) && (stepsize = 1.0)
        if (max_val - min_val)/stepsize < max_lines
            stepsize = max(round((max_val-min_val)/num_lines*10)/10,0.1)
            (stepsize == 0.3) && (stepsize = 0.2)
            (stepsize == 0.4) && (stepsize = 0.5)
            ((stepsize >= 0.6) && (stepsize <= 0.8)) && (stepsize = 0.5)
            (stepsize >= 0.9) && (stepsize = 1.0)
        end
        if (max_val - min_val)/stepsize < ceil(max_lines/2)
            stepsize = (max_val - min_val) / max_lines
        else
            min_val = ceil(min_val/stepsize) * stepsize
            max_val = floor(max_val/stepsize) * stepsize
        end
        if stepsize > 0
            levels = collect(min_val:stepsize:max_val)
        else
            levels = [min_val,max_val]
        end
    end
    ll = length(levels)
    levels = flipdim(levels[end:-max(floor(Int,ll/9),1):1],1)
    (length(levels) == 1) && (levels = levels .* ones(2))
    # upstream has commented-out normality message logic here
    return levels,err
end

"""
determine a grid size which will result in a quick computation
"""
function setgridsize(n,nmin,nmax,iscplx)
    npts = round(Int,min(max(300/n^(5/6),nmin),nmax))
    iscplx && (npts = max(nmin,floor(Int,3*npts/4)))
    isodd(npts) && (npts += 1)
    npts
end

"""
format a time interval (in seconds) into a short string
"""
function prettytime(ttime)
    if ttime < 180
        total_time = ceil(Int,ttime)
        str = @sprintf("%d seconds",total_time)
    elseif ttime < 10800
        total_time = ceil(Int,ttime/60)
        str = @sprintf("%d minutes",total_time)
    else
        total_time = round(Int,ttime/360)/10
        str = @sprintf("%.1f hours",total_time)
    end
    return total_time,str
end

isvalidax(ax) = false
function isvalidax{T<:Real}(ax::Vector{T})
    (length(ax) == 4) && (ax[1] < ax[2]) && (ax[3] < ax[4])
end

"""
construct a mesh for the imaginary axis to exploit symmetry

shift_axes(ax,npts) -> y,n_mirror
"""
function shift_axes(ax,npts)
    # note: linspace arg len::Float is ok
    if ax[4] > -ax[3]
        tpts = floor((ax[4]/(ax[4]-ax[3]))*npts)
        y = collect(linspace(ax[4]/(2*(tpts-1)+1),ax[4],tpts))
        n = count(x -> (x < -ax[3]),y)
        if abs(y[min(length(y),n+1)] + ax[3]) > 1e-15
            y = vcat(y[1:n],-ax[3],y[n+1:end])
        end
        num_mirror = n+1
    elseif ax[4] < -ax[3]
        bpts = floor((ax[3] / (ax[3] - ax[4]))*npts)
        y = collect(linspace(ax[3],ax[3]/(2*(bpts-1)+1),bpts))
        n = count(x -> (x < -ax[4]),y)
        if abs(y[max(1,n)] + ax[4]) > 1e-15
            y = vcat(y[1:n],-ax[4],y[n+1:end])
        end
        num_mirror = -(length(y)-(n+1)+1)
    else
        if isodd(npts)
            y = collect(linspace(0,ax[4],round(Int,(npts+1)/2)))
            num_mirror = ((npts+1) >> 1) -1
        else
            step = (ax[4] - ax[3]) / npts
            y = collect(linspace(step/2,ax[4],round(Int,npts/2)))
            num_mirror = npts >> 1
        end
    end
    return y,num_mirror
end

function get_step_size(n,ly,routine)
    if n < 100
        step = max(1,floor(Int,ly/8))
    else
        step = min(ly,max(1,floor(Int,4*ly/n)))
    end
    # upstream decreases by factor of 4 if fast implementation is missing
    return step
end

"""
Prompt user for key character, with optional explanatory details.
Returns index of selected option.
"""
function replqdlg(query,details=""; options=["Yes","No"])
    abbrevs = [s[1] for s in options]
    nopts = length(options)
    qstr = (isempty(details) ? "): " : ",?): ")
    while true
        print(query," (",join(options,","),qstr)
        str = readline()
        c = uppercase(str[1])
        if !isempty(details) && (c == '?')
            println(details)
            continue
        end
        for i in 1:nopts
            if c == abbrevs[i]
                return i
            end
        end
        println("Please respond with one of $(abbrevs).")
    end
    false
end

dummyqdlg(query,details=""; options=["Yes","No"]) = 1

"""
compute reasonable axis limits for displaying a vector of points & environs
"""
function vec2ax(dispvec)
    ca = [minimum(real(dispvec)), maximum(real(dispvec)),
          minimum(imag(dispvec)), maximum(imag(dispvec))]
    if ca[1]==ca[2]
        ca[1] -= 0.5
        ca[2] += 0.5
    end
    if ca[3]==ca[4]
        ca[3] -= 0.5
        ca[4] += 0.5
    end
    ext = max(ca[2]-ca[1], ca[4]-ca[3]) / 3
    ax = [ca[1]-ext,ca[2]+ext,ca[3]-ext,ca[4]+ext]
    ax[1],ax[2] = tidyaxes(ax[1],ax[2],0.02)
    ax[3],ax[4] = tidyaxes(ax[3],ax[4],0.02)
    return ax
end

"""
Compute triangular factor(s) for a non-square matrix.

As used to prepare for pseudospectra computation.
"""
function rect_fact(A)
    m,n = size(A)
    if m >= 2*n
        # QR form
        F = qrfact(A[n+1:end,:])
        S = vcat(A[1:n,:],F[:R])
        T = I
    else
        # QZ form
        tmp = eye(m,n)
        F = schurfact(A[end-n+1:end,:],tmp[end-n+1:end,:])
        S = vcat(A[1:m-n,:]*F[:Z],F[:S])
        T = vcat(F[:Z][1:m-n,:],F[:T])
    end
    return S,T
end

function expandlevels(ld::LevelDesc)
    if ld.isunif
        if ld.step != 0
            levels = ld.first:ld.step:ld.last
        else
            levels = ones(2)*ld.first
        end
    else
        levels = ld.full_levels
        (length(levels) == 1) && (levels = levels[1]*ones(2))
    end
    isa(levels,Range) && (levels = collect(levels))
    levels
end

"""
tidy up axis limits

    tidyaxes(vmin,vmax,tol) -> nmin,nmax

rounds [vmin,vmax] outward to have cleaner decimal expressions
"""
function tidyaxes(vmin::Number,vmax::Number,tol)
    span = vmax - vmin
    round_to = span * tol
    lround_to = log10(round_to)
    expon = floor(lround_to)
    mant = round(10^(lround_to-expon))
    round_to = mant*10^expon
    nmin = round_to*floor(vmin/round_to)
    nmax = round_to*ceil(vmax/round_to)

    return nmin,nmax
end
