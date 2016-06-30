# startup script

workspace()

include("linalg_helpers.jl")

include("CanonizableGraph.jl")

include( "SmallGraphs.jl" );

include( "SmallDiGraph.jl" );

include( "MultiGraph.jl" );

include("nauty_interface.jl")

include("SharedCode.jl");

include("GraphVectorSpace.jl");

include("GraphOperator.jl");

include("DisplayHelpers.jl");

include("OrdinaryGraphComplex.jl");

include("EdgeMarkedGraphComplex.jl");

include("RibbonGraphComplex.jl")
