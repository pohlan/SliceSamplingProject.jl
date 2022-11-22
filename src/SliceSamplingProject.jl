__precompile__(false)
module SliceSamplingProject

    export single_var

    using Distributions

    function single_var(f,                                # function proportional to density
                        x0::Float64;                      # current point
                        w::Float64=1.0,                   # starting width of the interval I
                        m::Int64=4,                       # for step_out, limits the size of I to m*w
                        p::Int64=4,                       # for doubling, limits ...
                        iters::Int64=1000,
                        method::String="step_out")
        x  = zeros(iters)
        xk = x0
        for i = 1:iters
            y  = rand(Uniform(0, f(xk)))
            if method == "step_out"
                I  = step_out(f, xk, y, w, m)
                xk = shrinkage(;f, x0=xk, y, I,w,doubling=false)
            elseif method == "doubling"
                I  = doubling(f, xk, y, w, p)
                xk = shrinkage(;f, x0=xk, y, I,w,doubling=true)
            else
                error("Specified method doesn't exist.")
            end
            x[i] = xk
        end
        return x
    end

    function shrinkage(;f, x0, y, I, w,
                        doubling=false) # whether or not the doubling procedure is used; if true, requires extra procedure for accept()
        Lbar, Rbar = I
        maxiters   = 1000 # shouldn't be necessary, should be guaranteed to terminate
        for i = 1:maxiters
            U  = rand(Uniform(0,1))
            x1 = Lbar + U * (Rbar-Lbar)
            if y<f(x1) && accept(;f,x0,x1,y,w,I,doubling)
                 return x1
            elseif x1 < x0
                Lbar = x1
            else
                Rbar = x1
            end
        end
        error("No acceptable x1 found.")
        return nothing
    end

    function accept(;f, x0, x1, y, w, I, doubling=false)
        if !doubling
            return true
        end
        Lhat, Rhat = I
        D    = false
        while Rhat-Lhat > 1.1*w
            M = (Lhat+Rhat)/2
            if (x0<M && x1>=M) || (x0>=M && x1<M)
                D = true
            end
            if x1<M
                Rhat = M
            else
                Lhat = M
            end
            if D && y>=f(Lhat) && y>= f(Rhat)
                return false
            end
        end
        return true
    end

    function step_out(f, x0, y, w, m)
        U = rand(Uniform(0,1))
        L = x0 - w*U
        R = L + w
        V = rand(Uniform(0, 1))
        J = floor(m*V)
        K = (m-1) - J
        while J>0 && y<f(L)
            L = L-w
            J = J-1
        end
        while K>0 && y<f(R)
            R = R+w
            K = K-1
        end
        return (L, R)
    end

    function doubling(f, x0, y, w, p)
        U = rand(Uniform(0,1))
        L = x0 - w*U
        R = L + w
        K = p
        while K > 0 && (y<f(L) || y<f(R))
            V = rand(Uniform(0,1))
            if V < 0.5
              L = L - (R-L)
            else
              R = R + (R-L)
            end
            K = K-1
        end
        return (L, R)
    end
end # module
