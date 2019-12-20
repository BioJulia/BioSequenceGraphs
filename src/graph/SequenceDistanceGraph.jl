###
### Node and Link types for SequenceDistanceGraph
###

const NodeID = Int64

include("SDGNode.jl")
include("DistanceGraphLink.jl")

const LinksT = Vector{Vector{DistanceGraphLink}}

###
### Graph types
###

"""
The SequenceDistanceGraph is a representation of a genome assembly.
Sequences are contained in nodes, and the distances are represented by links.

A singe node represents a sequence *and* its reverse complement.
Every node has a correlative ID starting from 1.
For every node X in the graph, the negative ID -X is mapped to the reverse
complement of X. This mapping is virtual: Only one node is stored in the graph.
This is because every node has an orientaton: Each node has a positive end (+),
and a negative end (-).
So when a node is accessed with (or traversed by entering) the positive end the
node yields the stored sequence.
Conversely, when a node is accessed with (or traversed by entering) the negative
end the node yelds the reverse complement of the stored sequence.
In this way the positive end can be thought of as the sequence start, and the
negative end can be thought of as the sequence end.

A single distance between two sequences is represented as a single link.
Every link connects two node ends and contains a distance (they take the form
`([+, -]n1, [+, -]n2, [+, -]dist)`).
A link connects two node ends, and so the order of the signed nodes in the links
does not change the link.
If the distance in a link is negative, this represents an overlap between two
sequences. These overlaps must be "perfect overlaps".

This SequenceDistanceGraph type is only intended to be interacted with directly
by developers and people who know what they are doing. Most tasks an end user
wants to do to manipulate or query a graph can be 
"""
struct SequenceDistanceGraph{S<:BioSequence}
    nodes::Vector{SDGNode{S}}
    links::LinksT
end

"Shorthand for SequenceDistanceGraph"
const SDG = SequenceDistanceGraph

include("SequenceGraphPath.jl")

"Construct an empty sequence distance graph."
function SequenceDistanceGraph{S}() where {S<:BioSequence}
    return SequenceDistanceGraph{S}(Vector{SDGNode{S}}(), LinksT())
end

##
## Internal / not-nessecerily-safe
##

"""
    check_node_id(sg::SequenceDistanceGraph, i::NodeID)

The method responsible for checking that a node id used as input to a method
that queries or edits a `SequenceDistanceGraph` is a sensible value.
"""
@inline function check_node_id(sg::SequenceDistanceGraph, i::NodeID)
    if 0 < abs(i) ≤ n_nodes(sg)
        return true
    end
    @error "Sequence graph has no node with ID of $i"
end


"""
Get a **reference** to the vector of nodes in a graph `sg`.

!!! warning
    It is a bad idea to edit this vector yourself unless you know what you are
    doing.
"""
@inline nodes(sg::SequenceDistanceGraph) = sg.nodes

"""
    node_unsafe(sg::SequenceDistanceGraph, n::NodeID)

Get a refernece specific node from a sequence distance graph `sg` using its
correlative node id `n`.

!!! note
    `node_unsafe` accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the node for node 5.

!!! warning
    This method returns a **reference** to a graph's underlying `SDGNode`.
    NOT a copy! Messing with this will screw up your graph. So this method is
    not recommended or exported - unless you really know what you're doing.

!!! warning
    This method is explicitly marked unsafe for a reason.
    It does *zero* checking of the value it is passed as the node id.
    It also makes use of the `@inbounds` macro.
    This makes it faster, but you need to be 100% sure that your code won't try
    to call this with a bad id. So it is not exported or recommended
    - unless you really know what you're doing.
"""
@inline node_unsafe(sg::SequenceDistanceGraph, n::NodeID) = @inbounds nodes(sg)[abs(n)]

"""
    node(sg::SequenceDistanceGraph, n::NodeID)

Get a **reference** to a specific node from a sequence distance graph `sg` using
its correlative node id `n`.

!!! note
    `node` accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the links for node 5.

!!! warning
    This method returns a **reference** to a graph's underlying `SDGNode`.
    NOT a copy! Messing with this will screw up your graph. So this method is
    not recommended or exported - unless you really know what you're doing.
"""
@inline function node(sg::SequenceDistanceGraph, n::NodeID)
    check_node_id(sg, n)
    return node_unsafe(sg, n)
end

"""
    sequence_unsafe(sg::SequenceDistanceGraph, n::NodeID)

Get the reference to a node's underlying sequence object.

!!! warning
    This method is unsafe, no checking of node id's occurs and you get a
    reference to the node's sequence object - not a copy, so doing transformation
    operations (reverse_complement, setindex, etc.) to it will probably screw up
    the graph!
"""
@inline sequence_unsafe(sg::SequenceDistanceGraph, n::NodeID) = sequence(node_unsafe(sg, n))

"""
Get a **reference** to the vector of vectors of links in a graph `sg`.

!!! warning
    It is a bad idea to edit this vector yourself unless you know what you are
    doing.
"""
@inline links(sg::SequenceDistanceGraph) = sg.links

"""
    links_unsafe(sg::SequenceDistanceGraph, n::NodeID)

Get a **reference** to a vector storing all the links of a node in a
`SequenceDistanceGraph`, the node is specified using its correlative node id `n`.

!!! note
    `links_unsafe` accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the links for node 5.

!!! warning
    This method returns a **reference** to an underlying links vector that the
    `SequenceDistanceGraph` owns - NOT a copy!
    Messing with this will screw up your graph. So this method is not
    recommended or exported - unless you really know what you're doing.

!!! warning
    This method is explicitly marked unsafe for a reason.
    It does *zero* checking of the value it is passed as the node id.
    It also makes use of the `@inbounds` macro.
    This makes it faster, but you need to be 100% sure that your code won't try
    to call this with a bad id. So it is not exported or recommended
    - unless you really know what you're doing.
"""
@inline links_unsafe(sg::SequenceDistanceGraph, n::NodeID) = @inbounds links(sg)[abs(n)]

"""
    links(sg::SequenceGraph, n::NodeID)

    Get a **reference** to a vector storing all the links of a node in a
    `SequenceDistanceGraph`, the node is specified using its correlative node id `n`.

!!! note
    `links` accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the links for node 5.

!!! warning
    This method returns a **reference** to an underlying links vector that the
    `SequenceDistanceGraph` owns - NOT a copy!
    Messing with this will screw up your graph. So this method is not
    recommended or exported - unless you really know what you're doing.
"""
@inline function links(sg::SequenceDistanceGraph, n::NodeID)
    check_node_id(sg, n)
    return links_unsafe(sg, n)
end

##
## Public / safe
##

"Get the name of the graph. Defaults to the symbol :sdg."
@inline name(sg::SequenceDistanceGraph) = :sdg

## Nodes and sequences

"Get the number of nodes in the sequence distance graph `sg`."
@inline n_nodes(sg::SequenceDistanceGraph) = length(nodes(sg))

"Iterate over every node ID in the sequence distance graph `sg`."
@inline each_node_id(sg::SequenceDistanceGraph) = eachindex(nodes(sg))

"""
    sequence(sg::SequenceDistanceGraph, n::NodeID)

Get the full sequence of a node in a sequence distance graph using its
correlative node id `n`.

!!! note
    `sequence` accepts a NodeID that can be positive or negative.
    Nodes represent stretches of sequence in a canonical orientation, if you ask
    for for the sequence of say the third node, the positive node id 3
    (which denotes traversing the third node in the forward direction),
    gives you the canonical sequence. If you use the negative ID -3
    (which denotes traversing the third node in the reverse direction), you will
    get the reverse complement of the node's canonical (forward) sequence.

!!! note
    It is safe to modify the returned sequence without screwing up your graph,
    yet thanks to BioSequences.jl's copy on write system for LongSequences, data
    copying will only occur if nessecery. You get the best of both worlds.
"""
function sequence(sg::SequenceDistanceGraph, n::NodeID)
    check_node_id(sg, n)
    seqref = sequence_unsafe(sg, n)
    outseq = typeof(seqref)(seqref, 1:lastindex(seqref))
    if n < 0
        reverse_complement!(outseq)
    end
    return outseq
end

## Graph Topology

"""
    find_link(sg::SequenceDistanceGraph, src::NodeID, dst::NodeID)

Find and return the link that exists between a source node and a destination
node, using their correlative node ids.

In this instance, the IDs also denote the "ends" of a node, a link connects.

Recall that every node has an orientaton: Each node has a positive end (+),
and a negative end (-).
So when a node is accessed with (or traversed by entering) the positive end the
node yields the stored sequence.
Conversely, when a node is accessed with (or traversed by entering) the negative
end the node yelds the reverse complement of the stored sequence.

Thus `find_link(sg, -5, 1)` means you want to find a link in the graph that
allows you to exit the (-) end of node 5 (meaning you just traversed it in the 
canonical + --> - orientation), and enter the (+) end of node 1 (meaning you will
also traverse node one in the canonical + --> - orientation).

By contrast `find_link(sg, -5, -1)` means you want to find a link in the graph
that allows you to exit the (-) end of node 5 - again having just traverse it in
the + --> - orientation, but want to enter node 1 through it's (-) end - and thus
traversing node 1 in the + <-- - orientation, yielding the reverse complement of
the canonical sequence of node 1.

If a link connecting the two node ends is present, you get it. If not, you get
`nothing`. So make sure to check the output.

!!! note
    This method is safe and public, but not reeeeaaallly intended for the end
    user, as they have to worry about links between "ends" of a node.
    In contrast to the rest of this framework, which rather has the user thinking
    in node centric terms.
    e.g. "I am at node 5 which means I have traversed it forwards or I am at
    node -5, which means I have traversed node 5 backwards. From here I can visit
    node 1, -3, and 2, so I may go through nodes 1 and 2 forwards, and through 3
    backwards. The methods `get_next_nodes` and `get_previous_nodes` are more
    user-friendly and follow this more node-centric mindset.
"""
function find_link(sg::SequenceDistanceGraph, src::NodeID, dst::NodeID)
    query = DistanceGraphLink(src, dst)
    for l in links(sg, src)
        if l == query
            return l
        end 
    end
    return nothing
end

"""
    forward_links(sg::SequenceDistanceGraph, n::NodeID)

Get a vector of the links that are ahread of you, as you traverse node `n`,
continuing forward in you present direction of travel (`n` can be a positive or
negative node id).

The node id can be positive or negative. For example if you use a positive ID,
such as 5, this means you are traversing node 5 in the canonical orientation,
and so you will get the links that allow you to leave node 5 in the canonical
direction. If you used a negative ID such as -5, that means you are traversing
node 5 in the non-canonical orientation, and so you will get the links that
allow you to leave node 5 in the non-canonical direction.
"""
function forward_links(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{DistanceGraphLink}()
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forwards_from(link, n)
            push!(r, link)
        end
    end
    return r
end

"""
    backward_links(sg::SequenceDistanceGraph, n::NodeID)

Get a vector of the links that are behind of you, as you traverse node `n`,
continuing forward in you present direction of travel (`n` can be a positive or
negative node id).

The node id can be positive or negative. For example if you use a positive ID,
such as 5, this means you are traversing node 5 in the canonical orientation,
and so you will get the links that would allow you to enter node 5 in the
canonical direction. If you used a negative ID such as -5, that means you are
traversing node 5 in the non-canonical orientation, and so you will get the links
that would allow you to enter node 5 in the non-canonical direction.
"""
backward_links(sg::SequenceDistanceGraph, n::NodeID) = forward_links(sg, -n)

"""
    get_next_nodes(sg::SequenceDistanceGraph, n::NodeID)

Get the nodes you may visit next as you exit node `n`, maintaining your current
direction of travel (`n` can be a positive or negative node ID).

!!! note
    The node id can be positive or negative. For example if you use a positive ID,
    such as 5, this means you are traversing node 5 in the canonical orientation,
    and so you will get the nodes you may visit as you leave node 5 in the canonical
    direction. If you use a negative ID, such as -5, this means you are
    traversing node 5 in the non-canonical orientation, and so you will get the
    nodes you may visit as you leave node 5 in the non-canonical direction.

!!! note
    The list of nodes returned are also signed to denote direction: Say
    you got `[1, -2, 3]` as a result of `get_next_nodes(sg, 5)`. That means as
    you leave node 5 after traversing it in the canonical direction, you may
    proceed to travel through nodes 1 and 3 in the canonical direction, and node
    2 in the non-canonical direction.
"""
function get_next_nodes(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{NodeID}()
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forwards_from(link, n)
            push!(r, destination(link))
        end
    end
    return r
end

"""
    get_previous_nodes(sg::SequenceDistanceGraph, n::NodeID)

Get the nodes you may have previously been on before you enterd node `n`,
maintaining your current direction of travel (`n` can be a positive or negative
node ID).

!!! note
    The node id can be positive or negative. For example if you use a positive ID,
    such as 5, this means you are traversing node 5 in the canonical orientation,
    and so you will get the nodes you may have visited prior to entering node 5
    in the canonical direction. If you use a negative ID, such as -5, this means
    you are traversing node 5 in the non-canonical orientation, and so you will
    get the nodes you may have visited prior to entering node 5 in the
    non-canonical direction.

!!! note
    The list of nodes returned are also signed to denote direction: Say
    you got `[1, -2, 3]` as a result of `get_previous_nodes(sg, 5)`. That means
    prior to entering node 5 in the canonical direction, you may have traveled
    through nodes 1 and 3 in the canonical direction, and node 2 in the
    non-canonical direction.
"""
function get_previous_nodes(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{NodeID}()
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_backwards_from(link, n)
            push!(r, -destination(link))
        end
    end
    return r
end

"""
    find_tip_nodes!(result::Set{NodeID}, sg::SequenceDistanceGraph, min_size::Integer)

Get a set of IDs of all the nodes in the graph `sg` that count as "tips".

Here a tip node is defined as a node that only has a single neighbouring node
at one of it's ends, and no neighbouring nodes out of one of its other end.

The node's sequence must also be larger than `min_size` base pairs in length.

!!! note
    This method modifys a `result` input. So if you want to find tips in a
    graph repeatedly, you can resuse a `Set`, saving you some allocation costs.
"""
function find_tip_nodes!(result::Set{NodeID}, sg::SequenceDistanceGraph, min_size::Integer)
    empty!(result)
    for n in each_node_id(sg)
        nd = node(sg, n)
        if is_deleted(nd) || length(nd) > min_size
            continue
        end
        fwl = forward_links(sg, n)
        bwl = backward_links(sg, n)
        if length(fwl) == 1 && length(bwl) == 0
            if length(backward_links(sg, destination(first(fwl)))) > 1
                push!(result, n)
            end
        end
        if length(fwl) == 0 && length(bwl) == 1
            if length(forward_links(sg, -destination(first(bwl)))) > 1
                push!(result, n)
            end
        end
        if isempty(fwl) && isempty(bwl)
            push!(result, n)
        end
    end
    return result
end

"""
    find_tip_nodes(sg::SequenceDistanceGraph, min_size::Integer)

Get a set of IDs of all the nodes in the graph `sg` that count as "tips".

Here a tip node is defined as a node that only has a single neighbouring node
at one of it's ends, and no neighbouring nodes out of one of its other end.

The node's sequence must also be larger than `min_size` base pairs in length.
"""
function find_tip_nodes(sg::SequenceDistanceGraph, min_size::Integer)
    return find_tip_nodes!(Set{NodeID}(), sg, min_size)
end



















###
### Graph editing operations
###

"""
    add_node!(sg::SequenceDistanceGraph{S}, n::SDGNode{S}) where {S<:BioSequence}

Add a node to a sequence distance graph.

Returns the node ID used to access the new node added in the graph.

!!! warning
    Currently, we don't enforce the sequence in the node is canonical here.
    We just trust that it is canonical.

!!! note
    Adding a node to the graph does just that. After adding the node it still
    will not be linked to any other nodes.
"""
function add_node!(sg::SequenceDistanceGraph{S}, n::SDGNode{S}) where {S<:BioSequence}
    newlen = length(push!(nodes(sg),n))
    push!(links(sg), Vector{DistanceGraphLink}())
    return newlen
end

"""
   add_node!(sg::SequenceDistanceGraph{S}, seq::BioSequence) where {S<:BioSequence}

Add a sequence to a sequence distance graph as a node.

Returns the node ID used to access the new node added in the graph.

Can accept any sequence type and will attempt to coerce the input sequence to the
type required by the graph.

!!! warning
    Currently, we don't enforce the sequence in the node is canonical here.
    We just trust that it is canonical.

!!! note
    Adding a node to the graph does just that. After adding the node it still
    will not be linked to any other nodes.
"""
function add_node!(sg::SequenceDistanceGraph{S}, seq::BioSequence) where {S<:BioSequence}
    return add_node!(sg, SDGNode{S}(convert(S, seq), false))
end

"""
    remove_node!(sg::SequenceDistanceGraph{S}, n::NodeID) where {S<:BioSequence}

Remove a node from a sequence distance graph.

!!! note
    This method can accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will end up deleting node 5.

!!! note
    Links involving this node will also be removed from the graph.
"""
function remove_node!(sg::SequenceDistanceGraph{S}, n::NodeID) where {S<:BioSequence}
    oldlinks = copy(links(sg, n))
    for oldlink in oldlinks
        remove_link!(sg, source(oldlink), destination(oldlink))
    end
    # TODO: This is a lazy solution to getting rid of the node.
    nodes(sg)[abs(n)] = empty_node(S)
end


"""
    add_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID, dist::Int)

Construct a link between two nodes in a sequence Graph.
"""
function add_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID, dist::Int)
    # Guard against someone adding links using ID's bigger than the current max NodeID in the graph.
    if abs(source) > length(links(sg))
        resize!(links(sg), abs(source))
    end
    if abs(dest) > length(links(sg))
        resize!(links(sg), abs(dest))
    end
    push!(links(sg, source), DistanceGraphLink(source, dest, dist))
    push!(links(sg, dest), DistanceGraphLink(dest, source, dist))
end

"""
    remove_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID)

Remove a link between two nodes in a SequenceDistanceGraph.
Returns a boolean indicating whether the removal was successful.
Reasons this function would not return `true` include that the link
didn't exist in the graph, and so could not be removed.
"""
function remove_link!(sg::SequenceDistanceGraph, src::NodeID, dst::NodeID)
    slinks = links(sg, src)
    slinkslen = length(slinks)
    filter!(!isequal(DistanceGraphLink(src, dst, 0)), slinks)
    dlinks = links(sg, dst)
    dlinkslen = length(dlinks)
    filter!(!isequal(DistanceGraphLink(dst, src, 0)), dlinks)
    return slinkslen != length(slinks) || dlinkslen != length(dlinks)
end

remove_link!(sg::SequenceDistanceGraph, lnk::DistanceGraphLink) = remove_link!(sg, source(lnk), destination(lnk))

"""
Removes all the links in the collection from and to a given nodeID.
"""
function disconnect_node!(sg::SequenceDistanceGraph, n::NodeID)
    for flink in forward_links(sg, n)
        remove_link!(sg, flink)
    end
    for rlink in backward_links(sg, n)
        remove_link!(sg, rlink)
    end
end

###
### Graph traversal
###








function write_to_gfa1(sg, filename)
    @info string("Saving graph to ", filename)
    fasta_filename = "$filename.fasta"
    gfa = open("$filename.gfa", "w")
    fasta = open(FASTA.Writer, fasta_filename)
    println(gfa, "H\tVN:Z:1.0")
    for nid in eachindex(nodes(sg))
        n = node(sg, nid)
        if n.deleted
            continue
        end
        println(gfa, "S\tseq", nid, "\t*\tLN:i:", length(n), "\tUR:Z:", fasta_filename)
        write(fasta, FASTA.Record(string("seq", nid), sequence(n)))
    end
    close(fasta)
    for ls in links(sg)
        for l in ls
            if source(l) <= destination(l)
                print(gfa, "L\t")
                if source(l) > 0
                    print(gfa, "seq", source(l), "\t-\t")
                else
                    print(gfa, "seq", -source(l), "\t+\t")
                end
                if destination(l) > 0
                    print(gfa, "seq", destination(l), "\t+\t")
                else
                    print(gfa, "seq", -destination(l), "\t-\t")
                end
                println(gfa, distance(l) < 0 ? -distance(l) : 0, "M")
            end
        end
    end
    close(gfa)
end

function add_nodes!(sg::SequenceDistanceGraph{S}, fa::FASTA.Reader) where {S<:BioSequence}
    rec = FASTA.Record()
    rcnodes = 0
    firstlen = n_nodes(sg)
    while !eof(fa)
        read!(fa, rec)
        if iscanonical(s)
            add_node!(sg, s)
        else
            add_node!(sg, reverse_complement!(s))
            rcnodes += 1
        end
    end
    @info string("Read ", n_nodes(sg) - firstlen, " nodes from file (", rcnodes, " canonised).")
    return sg
end



function find_all_unitigs!(unitigs::Vector{SequenceGraphPath{G}},
    sg::G, min_nodes::Integer) where {G<:SequenceDistanceGraph}
    empty!(unitigs)
    consumed = falses(n_nodes(sg))
    for n in each_node_id(sg)
        if consumed[n] || is_deleted(node(sg, n))
            continue
        end
        consumed[n] = true
        path = SequenceGraphPath(sg, [n])
        
        # Two passes, fw and bw, path is inverted twice, so still n is +
        for pass in 1:2
            fn = forward_links(sg, last(path))
            while length(fn) == 1
                dest = destination(first(fn))
                if (!consumed[abs(dest)]) && (length(backward_links(sg, dest)) == 1)
                    push!(path, dest)
                    consumed[abs(dest)] = true
                else
                    break
                end
                fn = forward_links(sg, last(path))
            end
            reverse!(path)
        end
        if n_nodes(path) >= min_nodes
            push!(unitigs, path)
        end
    end
    return unitigs
end

function find_all_unitigs(sg::G, min_nodes::Integer) where {G<:SequenceDistanceGraph}
    return find_all_unitigs!(Vector{SequenceGraphPath{G}}(), sg, min_nodes)
end

function collapse_all_unitigs!(unitigs::Vector{SequenceGraphPath{G}},
                               newnodes::Vector{NodeID},
                               sg::G,
                               min_nodes::Integer,
                               consume::Bool) where {G<:SequenceDistanceGraph}
    
    find_all_unitigs!(unitigs, sg, min_nodes)
    resize!(newnodes, length(unitigs))
    @inbounds for i in eachindex(unitigs)
        newnodes[i] = join_path!(unitigs[i], consume)
    end
end

function collapse_all_unitigs!(sg::SequenceDistanceGraph, min_nodes::Integer, consume::Bool)
    unitigs = Vector{SequenceGraphPath{typeof(sg)}}()
    newnodes = Vector{NodeID}()
    return collapse_all_unitigs!(unitigs, newnodes, sg, min_nodes, consume)
end

"""
    load_from_gfa1!(sg::SequenceDistanceGraph{S}, gfafile::AbstractString, fafile::AbstractString) where {S<:BioSequence}

Load a graph from a GFAv1 formatted file, and associated FASTA file of node sequences.

!!! note
    The GFA format permits storing sequences of the graph nodes in a seperate fasta
    file, instead of in the GFA file. This is so as the sequences of the graph
    nodes can be easily fed into other tools that typically accept FASTA files
    as input. Many assemblers also output a GFA + FASTA combo. 
    Therefore, this method asks for the filepath of a GFAv1 file, as well as a
    filepath to a FASTA formatted file. This method reads the node sequences from
    the FASTA file, before getting the links between nodes from the GFAv1
    file.   
"""
function load_from_gfa1!(sg::SequenceDistanceGraph{S},
                         gfafile::AbstractString,
                         fafile::AbstractString) where {S<:BioSequence}
    @info "Loading graph"
    @info string("Graph FASTA filename: ", fafile)
    # Load all the sequences from FASTA file. If they are not canonical, flip them,
    # and remember that they are flipped.
    @info string("Loading sequences from ", fafile)
    fardr = open(FASTA.Reader, fafile)
    add_nodes!(sg, fardr)
    @info string("Loading links from ", gfafile)
    gfas = open(gfafile, "r")
    for line in eachline(gfas)
        
    end
end

Base.summary(io::IO, sdg::SequenceDistanceGraph) = print(io, "Sequence distance graph (", n_nodes(sdg), " nodes)")
    