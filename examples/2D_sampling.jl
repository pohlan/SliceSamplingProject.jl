using SliceSamplingProject, LinearAlgebra, Distributions, PyPlot

μ1 = [0, 0]
μ2 = [10, 6]
μ3 = [2, -10]
#Σ = Matrix(3.0I, 3,3)
Σ = [0.5 0; 0 10.0]

f(x) = pdf(MvNormal(μ1, Σ), x) + pdf(MvNormal(μ2, 2*Σ),x) + pdf(MvNormal(μ3, 0.5*Σ),x)
#f(x) = 10 - (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
#f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 -7)^2
#f(x) = 0.01 + abs(sin(x[1]) * cos(x[2]) * exp(abs(1 - (sqrt(x[1]^2 + x[2]^2)/pi))))

x0 = [0.,0.]
x = slice_sampling_2D(f,x0,w=50.0,iters=10000)
# too small w requires higher # iterations

close("all")
h = hist2D(x[1,:],x[2,:],bins=1000)
