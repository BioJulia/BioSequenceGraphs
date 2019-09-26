var documenterSearchIndex = {"docs":
[{"location":"DeBruijnGraph/#de-Bruijn-Graph-type-1","page":"-","title":"de Bruijn Graph type","text":"","category":"section"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"A fundamental approach for de-novo gene assembly is to use make use of de Bruijn graphs. The graph is used to represent fragments of reads (mostly starting with kmers) as vertices and overlaps between these fragments as edges. DeBruijnGraph type is a special type of SequenceGraph. It is also made up of two fields:","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"nodes\nlinks","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"struct DeBruijnGraph\n    nodes::Vector{SequenceGraphNode}\n    links::Vector{Vector{SequenceGraphLink}}\nend","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"We initialize a DeBruijnGraph using the deBruijn_constructor function. This is mainly due to the fact that arbitrary links between two vertices are not allowed in the de Bruijn graph formalism. The constructor receives as input a list of kmers and generates the deBruijn_Graph where each kmer is a unique vertex and each overlap  of length k-1 is represented with an edge.","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"kmer_vector = generate_random_kmers(DNA,4,10)\ndbg = deBruijn_constructor(kmer_vector)\n\nDeBruijnGraph(SequenceGraphNode[SequenceGraphNode{Kmer{DNA,4}}(CGCC, true), SequenceGraphNode{Kmer{DNA,4}}(TCTG, true), SequenceGraphNode{Kmer{DNA,4}}(TGTG, true), SequenceGraphNode{Kmer{DNA,4}}(GAAG, true), SequenceGraphNode{Kmer{DNA,4}}(GGCA, true), SequenceGraphNode{Kmer{DNA,4}}(ACGA, true), SequenceGraphNode{Kmer{DNA,4}}(CGCT, true), SequenceGraphNode{Kmer{DNA,4}}(TCTC, true), SequenceGraphNode{Kmer{DNA,4}}(TACG, true), SequenceGraphNode{Kmer{DNA,4}}(GCAT, true)], Array{SequenceGraphLink,1}[[], [], [], [], [SequenceGraphLink(-5, 10, 1)], [], [], [], [SequenceGraphLink(-9, 6, 1)], []])","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"We can use the fastq readers to work on real data. The example below generates the set of unique kmers from 5-long reads. The kmer size is set to 15.","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"r = FASTQ.Reader(open(\"URnano_ecoli.fastq\", \"r\"))\nEcoli_reads = Set{BioSequence{DNAAlphabet{4}}}()\n\n## get first 5 reads\nfor i in 1:5\n    next_seq = iterate(r)\n    seq = sequence(next_seq[1])\n    next_seq = iterate(r,next_seq[2])\n    push!(Ecoli_reads,seq)\nend\nkmers = new_extract_canonical_kmers(Ecoli_reads,15)","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"The output is a set of kmers in their canonical form:","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"Set(Kmer{DNA,15}[AAAATATCTCGTTTT, TAAACCAGTCGCCGC, TTCGACATTACCCAG, CGCCCTGCCAGCAGT, TAATATTGTTCCATT, TGGTAATGGTCACAG, AAAAATTAAGCAGGA, ATATAAGTTATATCA, GCCCGATCTGTCTCC, GATTTCTCCGGGCCA  …  TGAGCGATTGCCTGA, CGGAGCAGCAGTGTC, AAAATCGTACATACC, GTCGCCTGATGCCTG, GTCAGCGAACCTTCC, CGCCGCTCACCGCCG, TGGATGAACGTTCAT, AATATGTCACAATTT, ATGCGATAGCAGGGG, GTAGAAAGCTCGTGG, CGATTGGTTTAAGAC])","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"Example data for checking node merging :","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"kmer_vector2 = Vector([DNAKmer{4}(\"ATTC\"),DNAKmer{4}(\"TTCG\"),DNAKmer{4}(\"TCGT\"),\n        DNAKmer{4}(\"AATC\"),DNAKmer{4}(\"AATG\"),DNAKmer{4}(\"CGTA\"),DNAKmer{4}(\"CGTC\")])","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"DeBruijnGraph for the above kmers:","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"DeBruijnGraph(Dict{Int64,SequenceGraphNode}(7=>SequenceGraphNode{Kmer{DNA,4}}(ACGA, true),4=>SequenceGraphNode{Kmer{DNA,4}}(CGTC, true),2=>SequenceGraphNode{Kmer{DNA,4}}(AATG, true),3=>SequenceGraphNode{Kmer{DNA,4}}(ATTC, true),5=>SequenceGraphNode{Kmer{DNA,4}}(CGTA, true),6=>SequenceGraphNode{Kmer{DNA,4}}(CGAA, true),1=>SequenceGraphNode{Kmer{DNA,4}}(AATC, true)), Dict(7=>[SequenceGraphLink(-7, 6, -3)],4=>[SequenceGraphLink(4, 7, -3)],2=>[],3=>[SequenceGraphLink(3, 1, -3), SequenceGraphLink(3, 2, -3)],5=>[SequenceGraphLink(5, 7, -3)],6=>[SequenceGraphLink(-6, -3, -3)],1=>[]), 4)","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"Then we apply node merging:","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"DeBruijnGraph(Dict{Int64,SequenceGraphNode}(7=>SequenceGraphNode{BioSequence{DNAAlphabet{4}}}(ACGAAT, true),4=>SequenceGraphNode{Kmer{DNA,4}}(CGTC, true),2=>SequenceGraphNode{Kmer{DNA,4}}(AATG, true),5=>SequenceGraphNode{Kmer{DNA,4}}(CGTA, true),1=>SequenceGraphNode{Kmer{DNA,4}}(AATC, true)), Dict(7=>[SequenceGraphLink(-7, 6, -3), SequenceGraphLink(-7, 1, -3), SequenceGraphLink(-7, 2, -3)],4=>[SequenceGraphLink(4, 7, -3)],2=>[],5=>[SequenceGraphLink(5, 7, -3)],1=>[]), 4)","category":"page"},{"location":"DeBruijnGraph/#","page":"-","title":"-","text":"So the nodes with nodeID 6 and 3 are collapsed into node 7 and we can see that both outgoing edges of 3 are given to node 7.","category":"page"},{"location":"types/graphs/#The-Sequence-Distance-Graph-(SDG)-1","page":"The Sequence Distance Graph (SDG)","title":"The Sequence Distance Graph (SDG)","text":"","category":"section"},{"location":"types/graphs/#","page":"The Sequence Distance Graph (SDG)","title":"The Sequence Distance Graph (SDG)","text":"SequenceDistanceGraph","category":"page"},{"location":"types/graphs/#GenomeGraphs.SequenceDistanceGraph","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.SequenceDistanceGraph","text":"The SequenceDistanceGraph is a representation of a genome assembly. Sequences are contained in nodes, and the distances are represented by links.\n\nA singe node represents a sequence and its reverse complement. Every node has a correlative ID starting from 1. For every node X in the graph, the negative ID -X is mapped to the reverse complement of X. This mapping is virtual: Only one node is stored in the graph. This is because every node has an orientaton: Each node has a positive end (+), and a negative end (-). So when a node is accessed with (or traversed by entering) the positive end the node yields the stored sequence. Conversely, when a node is accessed with (or traversed by entering) the negative end the node yelds the reverse complement of the stored sequence. In this way the positive end can be thought of as the sequence start, and the negative end can be thought of as the sequence end.\n\nA single distance between two sequences is represented as a single link. Every link connects two node ends and contains a distance (they take the form ([+, -]n1, [+, -]n2, [+, -]dist)). A link connects two node ends, and so the order of the signed nodes in the links does not change the link. If the distance in a link is negative, this represents an overlap between two sequences. These overlaps must be \"perfect overlaps\".\n\n\n\n\n\n","category":"type"},{"location":"types/graphs/#Querying-an-SDG-for-basic-properties-1","page":"The Sequence Distance Graph (SDG)","title":"Querying an SDG for basic properties","text":"","category":"section"},{"location":"types/graphs/#","page":"The Sequence Distance Graph (SDG)","title":"The Sequence Distance Graph (SDG)","text":"nodes\nn_nodes\neach_node_id\nlinks\nnode\nsequence","category":"page"},{"location":"types/graphs/#GenomeGraphs.nodes","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.nodes","text":"Get a reference to the vector of nodes in a graph sg.\n\nwarning: Warning\nIt is a bad idea to edit this vector yourself unless you know what you are doing. \n\n\n\n\n\n","category":"function"},{"location":"types/graphs/#GenomeGraphs.n_nodes","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.n_nodes","text":"Get the number of nodes in the sequence distance graph sg.\n\n\n\n\n\n","category":"function"},{"location":"types/graphs/#GenomeGraphs.each_node_id","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.each_node_id","text":"Iterate over every node ID in the sequence distance graph sg.\n\n\n\n\n\n","category":"function"},{"location":"types/graphs/#GenomeGraphs.links","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.links","text":"Get a reference to the vector of vectors of links in a graph sg.\n\nwarning: Warning\nIt is a bad idea to edit this vector yourself unless you know what you are doing.\n\n\n\n\n\nlinks(sg::SequenceGraph, n::NodeID)\n\nGet all of the links of a Node of a sequence distance graph using its  correlative node id n.\n\nnote: Note\nlinks accepts a NodeID that can be positive or negative. E.g. providing either 5 or -5 both mean node 5 in a graph, and so you will get the links for node 5.\n\n\n\n\n\n","category":"function"},{"location":"types/graphs/#GenomeGraphs.node","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.node","text":"node(sg::SequenceDistanceGraph, n::NodeID)\n\nGet a specific node from a sequence distance graph sg using its correlative node id n.\n\nnote: Note\nnode accepts a NodeID that can be positive or negative. E.g. providing either 5 or -5 both mean node 5 in a graph, and so you will get the links for node 5.\n\n\n\n\n\n","category":"function"},{"location":"types/graphs/#GenomeGraphs.sequence","page":"The Sequence Distance Graph (SDG)","title":"GenomeGraphs.sequence","text":"sequence(sg::SequenceDistanceGraph, n::NodeID)\n\nGet the full sequence of a node in a sequence distance graph using its correlative node id n.\n\n\n\n\n\n","category":"function"},{"location":"types/graphs/#Manually-editing-an-SDG-by-manipulating-nodes-and-links-1","page":"The Sequence Distance Graph (SDG)","title":"Manually editing an SDG by manipulating nodes and links","text":"","category":"section"},{"location":"types/graphs/#","page":"The Sequence Distance Graph (SDG)","title":"The Sequence Distance Graph (SDG)","text":"It is not recommended you do this if you are a high level user. However these small editing operations are required for developers. If you find yourself needing these methods, you will have to explicitly import them, as they are not exported from the module.","category":"page"},{"location":"types/graphs/#","page":"The Sequence Distance Graph (SDG)","title":"The Sequence Distance Graph (SDG)","text":"If you find yourself wanting to edit the graph manually, it's a good idea to ask the package authors listed in the Project.toml or .github/CODEOWNERS.","category":"page"},{"location":"types/graphs/#","page":"The Sequence Distance Graph (SDG)","title":"The Sequence Distance Graph (SDG)","text":"add_node!\nremove_node!\nadd_link!\nremove_link!\ndisconnect_node!","category":"page"},{"location":"man/guide/#Package-Guide-1","page":"Guide","title":"Package Guide","text":"","category":"section"},{"location":"man/guide/#Installation-1","page":"Guide","title":"Installation","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"GenomeGraphs is made available to install through BioJulia's package registry.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Julia by default only watches the \"General\" package registry, so before you start, you should add the BioJulia package registry.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Start a julia terminal, hit the ] key to enter pkg mode (you should see the prompt change from julia> to pkg>), then enter the following command:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"pkg> registry add https://github.com/BioJulia/BioJuliaRegistry.git","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"After you've added the registry, you can install GenomeGraphs from the julia REPL. Press ] to enter pkg mode again, and enter the following:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"pkg> add GenomeGraphs","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"man/guide/#Creating-a-WorkSpace-1","page":"Guide","title":"Creating a WorkSpace","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"You can create an empty genome graph workspace with the empty constructor:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"using GenomeGraphs\n\nws = WorkSpace()","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Obviously this workspace is quite useless on its own!","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"You need a genome graph and information to project onto it before you can do any exploration or analysis. There are a few ways to get your first graph:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Construct a de-Bruijn graph from raw sequencing reads.\nLoad a graph (such as produced from another assembler) from a GFAv1 file","category":"page"},{"location":"man/guide/#Constructing-a-de-Bruijn-graph-1","page":"Guide","title":"Constructing a de-Bruijn graph","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Let's see how to do option number 1, and construct a de-Bruijn graph from raw sequencing reads. This can be achieved with a few simple steps:","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Prepare the sequencing reads.","category":"page"},{"location":"man/guide/#Preparing-the-sequencing-reads-1","page":"Guide","title":"Preparing the sequencing reads","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Let's prepare the sequencing reads.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"WorkSpaces store sequencing reads in ReadDatastores, provided by the ReadDatastores.jl package.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"GenomeGraphs uses (and re-exports types and methods from) ReadDatastores. ReadDatastores is a standalone package in its own right (although it was built in the first place for GenomeGraphs).","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"If you want to use ReadDatastores as a standalone package in another Bio(Julia) based project (and we recommend you do - the data stores are more efficient than text files), you can find standalone docs for the package here. Some types and methods documented there, are repeated in the Library section of this manual here, for convenience.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Anyway, let's see how to build a paired end reads datastore!","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"using GenomeGraphs\n\nfwq = open(FASTQ.Reader, \"test/ecoli_pe_R1.fastq\")\nrvq = open(FASTQ.Reader, \"test/ecoli_pe_R2.fastq\")\n\nds = PairedReads(fwq, rvq, \"ecoli-test-paired\", \"my-ecoli\", 250, 300, 0, FwRv)","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Here \"ecoli-test-paired\" is provided as the base filename of the datastore, the datastore is given the name of \"my-ecoli\", this name will be used to identify it in the workspace later. The minimum length for the reads is set at 250 base pairs, and the maximum length is set to 300 base pairs. Reads that are too short will be discarded, reads that are too long are truncated.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"note: Note\nThe insert size of the paired reads to 0, since I'm not sure of it and right now the value is optional.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"I set the orientation of the paired reads to FwRv. This is the default, and means for every pair of reads, read 1 is oriented in the forward direction, and read 2 is oriented backwards (forwards on the opposite strand). This orientation distinguishes regular paired-end reads from other paired read types like Long Mate Pairs.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Now the datastore is created, it can be added to a workspace.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"ws = WorkSpace()\nadd_paired_reads!(ws, ds)","category":"page"},{"location":"man/guide/#Run-the-dbg-process-1","page":"Guide","title":"Run the dbg process","text":"","category":"section"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"GenomeGraphs comes with some very high-level methods in it's API, that we like to call processes. They perform some critical and common task as part of a larger workflow. Examples include constructing a de-Bruijn graph from sequencing reads, mapping reads to a graph, kmer counting and so on.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"Once a workspace has an attached read datastore, you can run the dbg process to produce a first de-Bruijn graph of the genome.","category":"page"},{"location":"man/guide/#","page":"Guide","title":"Guide","text":"dbg!(BigDNAMer{61}, 10, ws, \"my-ecoli\")","category":"page"},{"location":"types/workspace/#Workspaces-1","page":"Workspaces","title":"Workspaces","text":"","category":"section"},{"location":"#GenomeGraphs-1","page":"Home","title":"GenomeGraphs","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"A graph based genomics framework for the julia/BioJulia ecosystem.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: Latest release) (Image: MIT license)  (Image: Stable documentation) (Image: Pkg Status) (Image: Chat)","category":"page"},{"location":"#Introduction-1","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"GenomeGraphs is designed to do one thing - provide a framework that makes it simple for a human to work with genome graphs from scripts or interactively.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Graphs are the core representation used by genome assemblers to represent genome sequence models constructed from reads. At the time of writing it is fair to say that until recently, their use has been limited to the internals of genome assemblers, which are often treated as black boxes that output a series of flattened sequences in FASTA format.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The use of graphs has increased in recent years thanks to the GFA file format and developments in genome variation graphs and sequence to graph mappers.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"However, a lack of inter operation between graph-based tools, and limited tools for downstream graph-based analysis, contribute to a perceived complexity which maintains linear sequences and FASTA files as the typical unit of genomic sequence exchange.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Such flattening of graph representations within pipelines with multiple steps, that use different types of sequencing in an iterative fashion, produces ever-longer linear genome sequences through an information loss process.  As a result, genome assembly projects are prone to error propagation and difficult to reproduce and control.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"To address these problems for the BioJulia ecosystem, GenomeGraphs provides a flexible framework for building and integrating information over genome graphs.","category":"page"},{"location":"#Framework-overview-1","page":"Home","title":"Framework overview","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package implements a SequenceDistanceGraph type that represents sequences as nodes and the adjacency between sequences in links/edges. Rather than work directly with this graph data structure, you interact with a WorkSpace. A WorkSpace associates a SequenceDistanceGraph with raw sequencing reads,  sequence to graph mappings, and k-mer counts. The WorkSpace and the API  provides a working environment that enables you to project different kinds of  information over a graph, and navigate and analyse each node of a sequence graph.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Within the WorkSpace, you will find DataStore types permit random access to short, linked, and long read sequences stored on disk in BioSequences.jl's native bit encoding. Each datastore has an associated Mapper in the workspace that contains the output from mapping said reads onto the graph. KmerCounts allow you to compute k-mer coverage over the graph from sequencing data, enabling coverage analysis. Additional DistanceGraphs define alternative topologies over SequenceDistanceGraph nodes. They are typically used to represent longer range linkage information from various sequencing technologies.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Finally, a NodeView abstraction provides a proxy to a node, with methods to navigate a graph and access mapped data.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"This framework is intended to allow the user to expore genome graphs interactively and to create processing methods for assembly or downstream analyses.","category":"page"},{"location":"#TLDR;-Package-features-1","page":"Home","title":"TLDR; Package features","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Efficient on-disk (buffered) data stores for sequencing reads.\nA simple graph data-structure for representing assembled genomes.\nTODO: Transparent mapping of sequencing reads onto graphs.\nTODO: Kmer counts and coverage projection over genome graph nodes.\nWorkspaces binding a genome graph, mapped sequences, kmer counts, and annotation.\nDe-novo genome assembly utilities:\nde-Bruijn graph construction, with tip-clipping & bubble-popping.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"See the Guide section of the manual for a tutorial explaining how to get started using GenomeGraphs.","category":"page"}]
}
