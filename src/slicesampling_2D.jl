function slice_sampling_2D(f,                       # function proportional to density
                           x0;                      # current point
                           w=1.0,                   # starting width of the interval I
                           iters=1000)
    n  = length(x0)
    x  = zeros(n,iters)
    xk = x0
    for i = 1:iters
        y  = rand(Uniform(0, f(xk)))
        L = zeros(n); R = zeros(n)
        for i = 1:n
            # randomly position the hyperrectangle H
            U = rand(Uniform(0,1))
            L[i] = xk[i] - w*U
            R[i] = L[i] + w
        end
        maxiters = 1000
        for m = 1:maxiters
            # sample from H
            x1 = zeros(n)
            for i = 1:n
                U = rand(Uniform(0,1))
                x1[i] = L[i] + U * (R[i]-L[i])
            end
            if y < f(x1)
                xk = x1
                break
            end
            for i = 1:n
                if x1[i] < xk[i]
                    L[i] = x1[i]
                else
                    R[i] = x1[i]
                end
            end
        end
        x[:,i] = xk
    end
    return x
end
