module TrajectOpt

    using FLOWMath, SNOW, Plots, DifferentialEquations, StaticArrays, DelimitedFiles, Snopt
    DE = DifferentialEquations
    FM = FLOWMath
    SA = StaticArrays
    DF = DelimitedFiles
    SN = Snopt

    include("models.jl")
    include("dynamics.jl")
    include("optimize.jl")

    export Model
    export LowFidel 
    export HighFidel
    export BiWingTailSitter
    export polar_constructor
    export biwing_tailsitter_forces_constructor
    # export dynamicsEU
    # export simulateEU
    # export step!
    export dynamics2D!
    export simulate
    export plot_simulation
    export optimize_trim
    export optimize_trajectory
    export optimize_trajectory_by_segments

end
