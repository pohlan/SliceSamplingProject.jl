using SliceSamplingProject, PyPlot, Distributions
close("all") # close all open figures
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 28

fcts = Dict()

# Gaussian
σ1   = 2.0
fcts["Normal(3,2)"] = x -> pdf(Normal(3, 2), x)

# exponential
fcts["Exponential(1)"] = x -> pdf(Exponential(1.0), x)

fcts["Uniform(0,1)"] = x -> pdf(Uniform(0,1), x)

# Beta
fcts["Beta(3,10)"] = x -> pdf(Beta(3,10), x)

n  = length(fcts)
x0 = 0.5
figure(figsize=(28,9))
for (d, dist) in enumerate(keys(fcts))
    # stepping out procedure
    subplot(2,n,d)
    x = slice_sampling_1D(fcts[dist],x0;iters=10^5,w=1.0,method="step_out")
    hist(x,bins=40,color="grey",alpha=0.4,density=true,label="stepping out")
    xstart, xend = extrema(x)
    plot(xstart:0.01:xend,fcts[dist].(xstart:0.01:xend),color="black",lw=2,label="true pdf")
    title(dist)
    if d == 1 ylabel("f(x)") end
    if d == 2 legend() end
    # doubling
    subplot(2,n,d+n)
    x = slice_sampling_1D(fcts[dist],x0;iters=10^5,w=1.0,method="doubling")
    hist(x,bins=40,color="CornFlowerBlue",alpha=0.4,density=true,label="doubling")
    xstart, xend = extrema(x)
    plot(xstart:0.01:xend,fcts[dist].(xstart:0.01:xend),color="black",lw=2,label="true pdf")
    xlabel("x")
    if d == 1 ylabel("f(x)") end
    if d == 2 legend() end
end
savefig("simple_dists.jpeg")
