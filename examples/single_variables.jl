using SliceSamplingProject, PyPlot
close()

# simple Gaussian
σ    = 1.0
f(x) = 1/(σ*sqrt(2π)) * exp(-x^2/(2σ^2))
x0 = 1.0

x = single_var(f,x0;iters=50000)
hist(x,bins=50,color="grey",alpha=0.5,density=true);
plot(-4:0.1:4,f.(-4:0.1:4),"black",lw=3)
xlabel("x"); ylabel("f(x)")