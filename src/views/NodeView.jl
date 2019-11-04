

"The NodeView provides a read-only comfortable interface for graph traversal."
struct NodeView{G<:SequenceDistanceGraph}
    ws::WorkSpace
    id::NodeID
    graph::G
end

"""
The LinkView provides a read-only interface to links in the a SequenceDistanceGraph.

A LinkView can only be generated using the `next` and `prev` methods of the
`NodeView` type.
"""
struct LinkView{G<:SequenceDistanceGraph}
    node_view::NodeView{G}
    distance::Int32
end

"""
    Get a view of node `nodeid` in the base SequenceDistanceGraph of the
    WorkSpace `ws`.
"""
function node(ws::WorkSpace, nodeid::NodeID)
    return NodeView(ws, nodeid, ws.sdg)
end

node(lv::LinkView) = lv.node_view

"""
    Get the sequence of the underlying SequenceDistanceGraph node.
"""
sequence(nv::NodeView) = sequence(graph(nv), id(nv))

Base.length(nv::NodeView) = length(sequence(graph(nv), id(nv)))

id(nv::NodeView) = nv.id
graph(nv::NodeView) = nv.graph
workspace(nv::NodeView) = nv.ws

function next(nv::NodeView{G}) where {G}
    g = graph(nv)
    fwl = forward_links(g, id(nv))
    r = Vector{LinkView{G}}(undef, length(fwl))
    i = 1
    @inbounds for l in fwl
        r[i] = LinkView{G}(NodeView{G}(workspace(nv), destination(l), g), distance(l))
        i = i + 1
    end
    return r
end

function Base.summary(io::IO, nv::NodeView)
    print(io, "A view of a graph node (node: ", id(nv), ", graph: ", name(graph(nv)), "):")
end

function Base.show(io::IO, nv::NodeView)
    summary(io, nv)
    print(io, "\n  ")
    show(io, sequence(nv))
end

