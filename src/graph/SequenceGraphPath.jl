struct SequenceGraphPath{G<:SequenceDistanceGraph}
    sg::G
    nodes::Vector{NodeID}
end

SequenceGraphPath(sg::G) where {G<:SequenceDistanceGraph} = SequenceGraphPath{G}(sg, Vector{NodeID}())

@inline nodes(p::SequenceGraphPath) = p.nodes
@inline n_nodes(p::SequenceGraphPath) = length(nodes(p))
@inline Base.push!(p::SequenceGraphPath, n::NodeID) = push!(nodes(p), n)
@inline graph(p::SequenceGraphPath) = p.sg
@inline Base.last(p::SequenceGraphPath) = last(p.nodes)

function Base.reverse!(p::SequenceGraphPath)
    nds = nodes(p)
    for i in 1:div(lastindex(nds), 2)
        j = lastindex(nds) - i + 1
        @inbounds x = nds[i]
        @inbounds y = nds[j]
        @inbounds nds[i] = -y
        @inbounds nds[j] = -x
    end
    return p
end