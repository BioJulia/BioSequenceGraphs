__precompile__()

module GenomeGraphs

export
    ###
    ### Re-exports of GenomeGraphs framework sub-components.
    ###
    
    # Re-exports of ReadDatastores
    ReadDatastore,
    PairedReads,
    LongReads,
    LinkedReads,
    buffer,
    FwRv,
    
    # Re-exports of FASTX
    FASTA,
    FASTQ,
    
    # Re-exports of BioSequences
    DNAAlphabet,
    BioSequence,
    LongSequence,
    AbstractMer,
    Mer,
    BigMer,
    DNAMer,
    DNAKmer,
    BigDNAMer,
    BigDNAKmer,
    
    ###
    ### MerCounts
    ###
    #MerCounts, # Moved to KmerAnalysis.jl
    
    ### WorkSpace
    WorkSpace,
    add_paired_reads!,
    paired_reads,
    add_mer_counts!,
    mer_counts,
    
    ###
    ### Processes
    ###
    dbg,
    dbg!,
    remove_tips!

include("Graphs.jl")       # Submodule defining the key Graph type and basic methods.
include("GraphIndexes.jl") # Submodule defining types that allow indexing of a graph.

using BioSequences, FASTX, ReadDatastores, KmerAnalysis
import BioSequences.EveryMerIterator

include("workspace/WorkSpace.jl")
include("views/NodeView.jl")

include("processes/dbg.jl")
include("processes/remove_tips.jl")
end # module GenomeGraphs
