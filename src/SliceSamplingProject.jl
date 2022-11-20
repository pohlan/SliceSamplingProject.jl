module SliceSamplingProject

    using Distributions

    function single_var(f,x0)
        y  = rand(Uniform(0, f(x0)))
        I  = step_out(f, x0, y, w, m)
        x1 = shrinkage(;f, x0, y, I)
        return
    end

    function shrinkage(;f::Function,   # function proportional to density
                       x0::Float64,   # current point
                       y::Float64,    # vertical level defining the slice
                       I::Tuple{Float64, Float64})
        Lbar, Rbar = I
        for i = 1:maxiters
            U  = rand(Uniform(0,1))
            x1 = Lbar + U * (Rbar-Lbar)
            if y<f(x1) && accept(;f,x0,x1,y,w,L,R)
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

    # check again!!
    function accept(;f, x0, x1, y, w, L, R)
        Lhat = L
        Rhat = R
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
        return D
    end

    function step_out(f::Function,   # function proportional to density
                      x0::Float64,   # current point
                      y::Float64,    # vertical level defining the slice
                      w::Float64,    # estimate of the typical size of a slice
                      m::Int64)      # limiting the size of a slice to m*w
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

    function doubling(f::Function,   # function proportional to density
                      x0::Float64,   # current point
                      y::Float64,    # vertical level defining the slice
                      w::Float64,    # estimate of the typical size of a slice
                      p::Int64)      # limiting the size of a slize to 2^p*w
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
