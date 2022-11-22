using SliceSamplingProject, PyPlot, Distributions
close()

# Gaussian
σ1   = 1.0
#f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2))

# bimodal
σ2 = 3.0
# f(x) = 1/(σ1*sqrt(2π)) * exp(-x^2/(2σ1^2)) + 1/(σ2*sqrt(2π)) * exp(-(x-20)^2/(2σ2^2))

# exponential
# f(x) = pdf(Exponential(0.1), x)

f(x) = pdf(Uniform(0,1), x)

# Beta
# f(x) = pdf(Beta(1.0,10.0), x)

# run the slice sampling
x0 = 0.5
method = ["step_out", "doubling"][2]
x = single_var(f,x0;iters=100000,method)
hist(x,bins=50,color="grey",alpha=0.4,density=true);

# plotting
xstart,xend = extrema(x)
xs = range(xstart,xend,length=1000)
dx = xs[end]-xs[end-1]
k  = sum(f.(xs))*dx
plot(xs,f.(xs)./k,"black",lw=2.5)
xlabel("x"); ylabel("f(x)")
