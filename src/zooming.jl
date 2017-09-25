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
    zoomin!(ps_data,z,ax) -> Int

update `ps_data` to include a finer resolution portrait.

Returns -1 for invalidity, 0 for return to an existing zoom (requiring
a redraw), and 1 for establishment of new one (requiring a recalculation
and redraw).
"""
function zoomin!(ps_data::PSAStruct,z,ax)
    zr,zi = real(z[1]),imag(z[1])
    if ((zr < ax[1]) || (zr > ax[2]) || (zi < ax[3]) || (zi > ax[4]))
        return -1
    end
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    if length(z) == 1
        zr2,zi2 = zr,zi
    else
        zr2,zi2 = real(z[2]),imag(z[2])
    end
    zax = zoom.ax
    final_rect = [max(zax[1],min(zr,zr2)),min(zax[2],max(zr,zr2)),
                  max(zax[3],min(zi,zi2)),min(zax[4],max(zi,zi2))]
    aw = ax[2] - ax[1]
    ah = ax[4] - ax[3]
    # if rectangle is trivial, zoom in on first point, mag factor 2
    if (abs(final_rect[2]-final_rect[1])/aw < 0.01) ||
        (abs(final_rect[4]-final_rect[3])/ah < 0.01)
        # if at end of zoom chain, append
        if ps_data.zoom_pos == length(ps_data.zoom_list)
            w = zax[2]-zax[1]
            h = zax[4]-zax[3]
            # start with a copy of current portrait
            push!(ps_data.zoom_list,deepcopy(zoom))
            ps_data.zoom_pos += 1
            newzoom = ps_data.zoom_list[ps_data.zoom_pos]
            # work out new axes
            newax = [final_rect[1]-w/4,final_rect[1]+w/4,
                     final_rect[3]-h/4,final_rect[3]+h/4]
            copy!(newzoom.ax,newax)
            final_rect = newzoom.ax # note: an alias

            # for real matrix, if box is almost up-down symmetric,
            # enforce symmetry
            if ps_data.ps_dict[:Aisreal]
                span = final_rect[4] - final_rect[3]
                the_diff = (final_rect[4] + final_rect[3])/span
                if abs(the_diff) < 0.15
                    final_rect[3] = -max(abs.(final_rect[3:4]))
                    final_rect[4] = max(abs.(final_rect[3:4]))
                end
            end
            ax1 = newzoom.ax
            ax1[1],ax1[2] = tidyaxes(ax1[1],ax1[2],0.05)
            ax1[3],ax1[4] = tidyaxes(ax1[3],ax1[4],0.05)
            retval = 1
        else # there is already a zoom to go to
            ps_data.zoom_pos += 1
            retval = 0
        end
    else # nontrivial rectangle branch
        # start with a copy of current zoom
        push!(ps_data.zoom_list,deepcopy(zoom))
        ps_data.zoom_pos += 1
        # drop any obsolete zooms
        resize!(ps_data.zoom_list,ps_data.zoom_pos)
        newzoom = ps_data.zoom_list[ps_data.zoom_pos]
        # for real matrix, if box is almost up-down symmetric, enforce symmetry
        if ps_data.ps_dict[:Aisreal]
            span = final_rect[4] - final_rect[3]
            the_diff = (final_rect[4] + final_rect[3])/span
            if abs(the_diff) < 0.15
                final_rect[3] = -max(abs.(final_rect[3:4]))
                final_rect[4] = max(abs.(final_rect[3:4]))
            end
        end
        copy!(newzoom.ax,final_rect)
        ax1 = newzoom.ax
        ax1[1],ax1[2] = tidyaxes(ax1[1],ax1[2],0.02)
        ax1[3],ax1[4] = tidyaxes(ax1[3],ax1[4],0.02)
        newzoom.computed = false
        retval = 1
    end
    retval
end

"""
    zoomout!(ps_data,z,ax) -> Int

update `ps_data` to include a coarser resolution portrait.

Returns -1 for invalidity, 0 for return to an existing zoom (requiring
a redraw), and 1 for establishment of new one (requiring a recalculation
and redraw).
"""
function zoomout!(ps_data::PSAStruct,z,ax; include_fov=false)
    zr,zi = real(z[1]),imag(z[1])
    if (zr < ax[1]) || (zr > ax[2]) || (zi < ax[3]) || (zi > ax[4])
        return -1
    end
    if ps_data.zoom_pos > 1
        # restore previous zoom level
        ps_data.zoom_pos -= 1
        retval = 0
    else # no previous zoom entry to use; zoom out by factor of 2
        newzoom = deepcopy(ps_data.zoom_list[1])
        newzoom.computed = false
        unshift!(ps_data.zoom_list, newzoom)
        zax = newzoom.ax
        ps_dict = ps_data.ps_dict
        if include_fov && !isempty(get(ps_dict,:fov,zeros(0))) &&
              get(ps_dict,:show_fov,false)
            mnr,mxr = extrema(real(ps_dict[:fov]))
            mni,mxi = extrema(imag(ps_dict[:fov]))
            zoom2fov = (mnr <= zax[1]) || (mxr >= zax[2]) || (mni <= zax[3]) ||
                (mxi >= zax[4])
        else
            zoom2fov = false
        end
        if zoom2fov
            w = mxr - mnr
            h = mxi - mni
            ps_data.zoom_pos = 1
            zax[1] = mnr - w/4
            zax[2] = mxr + w/4
            zax[3] = mni - h/4
            zax[4] = mxi + h/4
        else
            # grow around selected point
            w = zax[2] - zax[1]
            h = zax[4] - zax[3]
            ps_data.zoom_pos = 1
            zax[1] = zr - w
            zax[2] = zr + w
            zax[3] = zi - h
            zax[4] = zi + h
        end
        zax[1],zax[2] = tidyaxes(zax[1],zax[2],0.05)
        zax[3],zax[4] = tidyaxes(zax[3],zax[4],0.05)
        retval = 1
    end
    retval
end
