using SliceSamplingProject, PyPlot, Distributions, StatsBase, LinearAlgebra, Infiltrator
close()

n = 120

# Gaussian
σ1   = 2.0
#f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2))

# bimodal
σ2 = 12.0
f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2)) + 1/(σ2*sqrt(2π)) * exp(-(x-60)^2/(2σ2^2))
xstart,xend = (-15., 100.)
xs = range(xstart,xend,length=n)
dx = xs[2]-xs[1]
xedges = range(xstart-0.5dx, xend+0.5dx, length=n+1)
k  = sum(f.(xs))*dx
true_pdf = f.(xs)./k


# exponential
# f(x) = pdf(Exponential(0.1), x)

# f(x) = pdf(Uniform(0,1), x)

# Beta
# f(x) = pdf(Beta(1.0,10.0), x)

# run the slice sampling
x0  = 0.5
ws = 0.1:0.2:10.0
its = floor.(Int,10 .^(2:0.2:6))
methods = ["step_out", "doubling"]
# err = zeros(length(its),length(methods))
err = zeros(length(ws),length(methods))
hs  = zeros(n,length(methods),length(ws))
true_pdf  = zeros(n,length(methods),length(ws))
for (m, method) in enumerate(methods)
    for (i,w) in enumerate(ws)
        # method = "step_out"; iters = 1000000
        x = single_var(f,x0;iters=10^4,w,method)
        figure(1)
        h_hist = hist(x,bins=n,color="grey",alpha=0.4,density=true);
        # h = fit(Histogram, x, xedges)
        # h = normalize(h, mode=:density)
        # err[i,m] = norm(h.weights - f.(xs)./k)
        err[i,m] = norm(h_hist[1] - f.(xs)./k,1)
        # hs[:,m,i] = h.weights
        hs[:,m,i] = h_hist[2][2:end]
    end
end
figure(2)
plot(ws,err[:,1],label="step out")
plot(ws,err[:,2],label="doubling")
title("$n bins")
legend()
xscale("log")
yscale("log")

# for b = 80:40:140
#         figure(1); h2 = hist(x,bins=b,density=true)
#         figure(2); plot(h2[2][2:end],h2[1],label="$b")
# end
# legend()

# x_p = sample(xs,Weights(f.(xs)./k),n; ordered=true)


# Metropolis
ws  = 1.0:5:30
err = zeros(length(ws))
hs  = zeros(n,length(ws))
acc_ratios = zeros(length(ws))
#     for (i,iters) in enumerate(its)
iters = 10^4
for (i,MH_width) in enumerate(ws)
    x, acc = metropolis(x0, x->log(f(x)) ; nsims=iters, MH_width)
    figure(1)
    h_hist = hist(x,bins=n,color="grey",alpha=0.4,density=true);
    err[i] = norm(h_hist[1] - f.(xs)./k,1)
    hs[:,i] = h_hist[2][2:end]
    acc_ratios[i] = acc
end

figure(13)
plot(ws,err)
title("Metropolis, err")
yscale("log")
figure(14)
plot(ws,acc_ratios)
title("Metropolis, acc_ratios")
