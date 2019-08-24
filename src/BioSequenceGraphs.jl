__precompile__()

module BioSequenceGraphs

export
    ###
    ### Sequence Distance Graph
    ###
    SequenceDistanceGraph,
    
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
    
    ###
    ### MerFreq
    ###
    collapse_into_freqs!,
    collapse_into_freqs,
    collapse_into_freqs_sorted!,
    collapse_into_freqs_sorted,
    merge_into!,
    merge_into_sorted!,
    
    
    
    ### WorkSpace
    WorkSpace,
    add_paired_reads!

using BioSequences, FASTX, ReadDatastores

include("mertools/MerFreq.jl")

include("graph/SequenceDistanceGraph.jl")
include("graph/graph_building.jl")
include("workspace/WorkSpace.jl")
end # module BioSequenceGraphs
