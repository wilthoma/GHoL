# startup script

workspace()

include("Settings.jl")

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

include("OrdinaryGraphComplexWrapper.jl");

include("HairyGraphComplex.jl")

include("EdgeMarkedGraphComplex.jl");

include("RibbonGraphComplex.jl")
