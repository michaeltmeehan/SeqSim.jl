"""
    TreeSequenceSimulation

Layered output for tree-driven sequence simulation.

`node_sequences` is the primary all-node simulation result. For TreeSim
adaptors, it is indexed directly by `TreeSim.Tree` node id, matching TreeSim's
canonical node-id-as-vector-index contract. `tip_alignment` is a derived
terminal projection of the same simulation. `root` is the tree root id used for
propagation, and `tips` is the vector of tree tip ids in the same order as
`tip_alignment`.
"""
struct TreeSequenceSimulation
    node_sequences::Vector{Sequence}
    tip_alignment::Alignment
    root::Int
    tips::Vector{Int}
end

"""
    simulate_on_tree(...)
    simulate_tree_sequences(...)

Simulate sequences over a tree-like object.

Concrete tree package adaptors define methods for this explicit bridge hook.
SeqSim's core intentionally provides only the API boundary and result type.
`simulate_on_tree` is retained as a compatibility alias where concrete adaptors
provide it; prefer `simulate_tree_sequences` for new all-node workflows.
"""
function simulate_on_tree end

"""
    simulate_tree_sequences(...)

Simulate sequences over a tree-like object and return a
`TreeSequenceSimulation`.

This is the preferred public name when callers need the all-node result. Package
extensions provide concrete methods for supported tree packages.
"""
function simulate_tree_sequences end

"""
    simulate_alignment(...)

Simulate sequences over a tree-like object and return the tip-only alignment.

This is a thin convenience surface for workflows that only need terminal
sequences. Package extensions provide concrete methods for supported tree
packages.
"""
function simulate_alignment end
