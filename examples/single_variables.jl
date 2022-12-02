using SliceSamplingProject, PyPlot, Distributions, StatsBase, LinearAlgebra
close("all") # close all open figures

n = 150

# Gaussian
σ1   = 2.0
#f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2))

# exponential
# f(x) = pdf(Exponential(1.0), x)
# xstart,xend = (0.,10.)
# k = 1.

# bimodal
σ2 = 12.0
f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2)) + 1/(σ2*sqrt(2π)) * exp(-(x-60)^2/(2σ2^2))
xstart,xend = (-15., 100.)
k = 2.

xs = range(xstart,xend,length=n)
dx = xs[2]-xs[1]
xedges = range(xstart, xend + dx, length=n+1)
true_pdf = f.(xs) ./k


# f(x) = pdf(Uniform(0,1), x)

# Beta
# f(x) = pdf(Beta(1.0,10.0), x)

# run the slice sampling
x0  = 0.5

ws = 1.0:1.:12.0
its = floor.(Int,10 .^[4,5,6,7])
err = Dict{String, Any}("step_out" => zeros(length(ws),length(its)),
                        "doubling" => zeros(length(ws),length(its)))
hs  = Dict{String, Any}("step_out" => zeros(length(ws),length(its),n),
                        "doubling" => zeros(length(ws),length(its),n))

# Slice sampling
for (m, method) in enumerate(keys(err))
    for (j,iters) in enumerate(its)
        println("Slice sampling for $iters iterations..")
        for (i,w) in enumerate(ws)
            x = single_var(f,x0;iters,w,method)
            h = fit(Histogram, x, xedges)
            h_density = h.weights ./ (sum(h.weights)*dx)            # specific for exponential function!
            err[method][i,j] = norm(h_density - true_pdf,1)
            hs[method][i,j,:] = h_density
        end
    end
end

# Metropolis
MHs  = ws
err["Metropolis"] = zeros(length(ws),length(its))
hs["Metropolis"]  = zeros(length(ws),length(its),n)
acc_ratios = zeros(length(ws))

for (j,iters) in enumerate(its)
    println("Metropolis sampling for $iters iterations..")
    for (i,MH_width) in enumerate(MHs)
        x, acc = metropolis(x0, x->log(f(x)) ; nsims=iters, MH_width)
        h = fit(Histogram, x, xedges)
        h_density = h.weights ./ (sum(h.weights)*dx)
        err["Metropolis"][i,j] = norm(h_density - true_pdf,1)
        hs["Metropolis"][i,j,:] = h_density
        acc_ratios[i] = acc
    end
end

figure(1)
for (m,method) in enumerate(keys(err))
    subplot(1,3,m)
    pcolor(err[method][2:end,:]); title(method); colorbar(); clim([0, maximum([maximum(err[m]) for m in keys(err)]) ])
    xlabel("# iterations"); ylabel("w")
end
