
"""
    get_mode(A, z, pseudo, verbosity)

Compute the (pseudo)mode of `A` associated with `z`.
"""
function get_mode(A, z, pseudo, verbosity)
    if pseudo
        n = size(A, 2)
        q = randn(n) + im*randn(n)
        normalize!(q)
        niter = 20

        SS, q = psmode_inv_lanczos(A, q, z, 1e-15, niter)   
    else
        SS, q = oneeigcond(A,z, verbosity)
    end
    return SS, q
end

"""
    get_mtxpowersnorm(A, nmax)

Compute ``|A^k|`` for k up to `nmax`.
"""
function get_mtxpowersnorm(A, nmax)
    transient = Vector{Float64}(undef, nmax+1)
    transient[1] = 1.0
    (min_tp, max_tp) = (Inf, -Inf)
    mat = copy(A)
    i = 2
    transient[i] = norm(mat)
    stop = false
    while !stop && i < nmax + 1
        mat = A * mat
        i += 1
        transient[i] = norm(mat)
        min_tp = min(min_tp, transient[i])
        max_tp = max(max_tp, transient[i])
        if (max_tp > 1e130) || ((min_tp < 1e-130) && (min_tp != 0.0))
            @warn("stopping: bounds going out of range")
            stop = true
        end
        if min_tp == 0.0
            @warn("stopping: exactly zero matrix")
            stop = true
        end
    end
    return 0:(i-1), transient[1:i]
end

"""
    get_mtxexpnorm(A, dt, nmax)

Compute ``|e^{A*k*dt}|`` for k up to `nmax`.
"""
function get_mtxexpnorm(A, dt, nmax)
    eAdt = exp(dt*A)
    ns, transient = get_mtxpowersnorm(eAdt, nmax)
    return dt*ns, transient
end