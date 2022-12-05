using SliceSamplingProject, PyPlot, Distributions, StatsBase, LinearAlgebra, LaTeXStrings
close("all") # close all open figures
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 30

n = 200

# bimodal distribution
σ1 = 2.0
σ2 = 12.0
f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2)) + 1/(σ2*sqrt(2π)) * exp(-(x-60)^2/(2σ2^2))
xstart,xend = (-15., 100.)
k = 2.

xs = range(xstart,xend,length=n)
dx = xs[2]-xs[1]
xedges = range(xstart, xend + dx, length=n+1)
true_pdf = f.(xs) ./k
figure(11,figsize=(12,7))
plot(xs,true_pdf,lw=2,"black")
xlabel("x"); ylabel("f(x)")
savefig("true_pdf.jpeg")

# run the slice sampling
x0  = 0.5
ws = 1.0:1.:45
its = floor.(Int,10 .^[4,4.5,5,5.5,6])
err = Dict{String, Any}("step_out" => zeros(length(ws),length(its)),
                        "doubling" => zeros(length(ws),length(its)))
hs  = Dict{String, Any}("step_out" => zeros(length(ws),length(its),n),
                        "doubling" => zeros(length(ws),length(its),n))

# Slice sampling
for (m, method) in enumerate(keys(err))
    for (j,iters) in enumerate(its)
        println("Slice sampling for $iters iterations..")
        for (i,w) in enumerate(ws)
            x = slice_sampling_1D(f,x0;iters,w,method)
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

# compare error of all
figure(1, figsize=(18,12))
for (m,method) in enumerate(keys(err))
    subplot(1,3,m)
    pcolor(log.(10,its),Array(ws),err[method][1:end,:]); title(method)
    xlabel(L"$log_{10}$ (iterations)")
    if m == 1
        ylabel(L"w")
    end
    if m == 3
        c = colorbar(); clim([0, maximum([maximum(err[m]) for m in keys(err)]) ])
        c.set_label("error")
    end
end
savefig("error_comparison.jpeg")

# some exapmles
figure(2, figsize=(16,10))
subplot(1,2,1)
for (m,method) in enumerate(keys(err))
    plot(xs,hs[method][15,3,:],lw=2.5,zorder=m,label=method)
end
plot(xs,true_pdf,"k",lw=2.5, label="true pdf")
xlabel("x"); ylabel("f(x)")
legend()
title(L"w = 15.0, $10^5$ iterations")

subplot(1,2,2)
for (m,method) in enumerate(keys(err))
    plot(xs,hs[method][15,5,:],lw=2.5,zorder=m,label=method)
end
plot(xs,true_pdf,"k",lw=2.5, label="true pdf")
xlabel("x")
title(L"w = 15.0, $10^6$ iterations")
savefig("pdf_comparison.jpeg")
