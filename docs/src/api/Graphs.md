# API: Graph submodule

## Types

```@docs
SequenceDistanceGraph
SDG
SDGNode
DistanceGraphLink
```

## Public / Safe methods

### Graph nodes and sequences

```@docs
name
n_nodes
sequence
```

### Graph Topology

```@docs
source
destination
distance
is_forwards_from
is_backwards_from
find_link
forward_links
backward_links
get_next_nodes
get_previous_nodes
find_tip_nodes
find_tip_nodes!
find_all_unitigs
find_all_unitigs!
```

### Graph IO

```@docs
write_to_gfa1
load_from_gfa1!
```

## Internal / Unsafe methods

### Graph nodes and sequences

```@docs
Graphs.empty_seq
Graphs.is_deleted
Graphs.sequence_unsafe
Graphs.check_node_id
Graphs.nodes
Graphs.node_unsafe
Graphs.node
```

### Graph topology

```@docs
Graphs.links
Graphs.links_unsafe
```

### Graph editing and manipulation

```@docs
Graphs.add_node!
Graphs.remove_node!
Graphs.add_link!
Graphs.remove_link!
Graphs.disconnect_node!
Graphs.collapse_all_unitigs!
```
