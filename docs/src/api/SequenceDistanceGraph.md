# API: Graph submodule

## Types

```@docs
SequenceDistanceGraph
SDG
DistanceGraphLink
SDGNode
```

### Public / Safe methods

#### Graph Nodes and sequences

```@docs
name
n_nodes
sequence
```

#### Graph Topology

```@docs
find_link
forward_links
backward_links
get_next_nodes
get_previous_nodes
find_tip_nodes
find_tip_nodes!
```

#### Graph IO

```@docs
write_to_gfa1
load_from_gfa1!
```

### Internal / Unsafe methods

#### Accessors and property queries
```@docs
GenomeGraphs.check_node_id
GenomeGraphs.nodes
GenomeGraphs.node_unsafe
GenomeGraphs.sequence_unsafe
GenomeGraphs.links
GenomeGraphs.links_unsafe
```