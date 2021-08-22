
_basic_psa_opts(zoom,ps_dict) = Dict{Symbol,Any}(
    #:levels=>expandlevels(zoom.levels),
    :recompute_levels=>zoom.autolev,
    :proj_lev=>zoom.proj_lev,
    :scale_equal=>zoom.scale_equal,
    :real_matrix=>ps_dict[:Aisreal],
    :verbosity=>0)