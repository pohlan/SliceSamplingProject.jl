__precompile__(false)
module SliceSamplingProject

    export slice_sampling_1D, slice_sampling_2D, metropolis

    using Distributions

    include("slicesampling_1D.jl")
    include("slicesampling_2D.jl")
    include("Metropolis.jl")

end # module
