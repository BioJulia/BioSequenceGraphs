__precompile__()

module BioSequenceGraphs

export
    ### Sequence Distance Graph
    SequenceDistanceGraph,
    DistanceGraphLink,
    SDGNode,
    # Basic queries and properties
    nodes,
    n_nodes,
    each_node_id,
    node,
    links,
    sequence,
    # Graph traversal
    get_next_nodes,
    get_previous_nodes,
    get_all_unitigs,
    build_unitigs_from_kmerlist2,
    delete_tips

using BioSequences

include("graph/SequenceDistanceGraph.jl")
include("graph/graph_building.jl")
end # module BioSequenceGraphs
