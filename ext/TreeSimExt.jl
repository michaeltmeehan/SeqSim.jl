module TreeSimExt

using Random
using SeqSim
using TreeSim

function _tree_taxon(tree::TreeSim.Tree, node::Int)
    label = tree.label[node]
    return label == 0 ? "node_$node" : label
end

function _validate_sequence_tree(tree::TreeSim.Tree)
    length(tree) > 0 || throw(ArgumentError("Cannot simulate sequences on an empty tree."))
    try
        TreeSim.validate_tree(tree)
    catch err
        throw(ArgumentError("TreeSim.Tree must be a validated canonical rooted tree: $(sprint(showerror, err))"))
    end
    return nothing
end

"""
    SeqSim.simulate_tree_sequences(rng, tree::TreeSim.Tree, site_model::SeqSim.SiteModel)
    SeqSim.simulate_tree_sequences(tree::TreeSim.Tree, site_model::SeqSim.SiteModel)

Simulate DNA sequences along a canonical `TreeSim.Tree`.

The root sequence is drawn from the substitution model frequencies in
`site_model`, then descendant node sequences are propagated in preorder using
`TreeSim.branch_length(tree, node)` as the nonnegative evolutionary length of
the incoming edge to each non-root node. Node times and labels are preserved as
sequence metadata; branch lengths and rooted structure are the simulation
inputs.

The returned `TreeSequenceSimulation` contains the primary all-node result in
`node_sequences`, indexed directly by TreeSim node id, plus the derived
tip-only projection in `tip_alignment`. Nonzero `tree.label[node]` values are
used as `Sequence.taxon`; nodes with label `0` use the deterministic fallback
`"node_\$(node)"`.
"""
function SeqSim.simulate_tree_sequences(
    rng::AbstractRNG,
    tree::TreeSim.Tree,
    site_model::SeqSim.SiteModel,
)::SeqSim.TreeSequenceSimulation
    _validate_sequence_tree(tree)

    root = TreeSim.root(tree)
    tip_ids = TreeSim.tips(tree)
    node_sequences = Vector{SeqSim.Sequence}(undef, length(tree))
    encoded_sequences = Vector{Vector{UInt8}}(undef, length(tree))

    root_sequence = SeqSim.rand_seq(rng, site_model; taxon=_tree_taxon(tree, root), time=tree.time[root])
    node_sequences[root] = root_sequence
    encoded_sequences[root] = SeqSim.encode(root_sequence.value)

    propagator = SeqSim.SequencePropagator(site_model)
    for node in Iterators.drop(TreeSim.preorder(tree, root), 1)
        parent = TreeSim.parent(tree, node)
        parent == 0 && throw(ArgumentError("Non-root node $node has no parent."))

        incoming_branch_length = TreeSim.branch_length(tree, node)
        encoded = propagator(rng, encoded_sequences[parent], incoming_branch_length)
        sequence = SeqSim.Sequence(encoded; taxon=_tree_taxon(tree, node), time=tree.time[node])
        encoded_sequences[node] = encoded
        node_sequences[node] = sequence
    end

    tip_alignment = [node_sequences[node] for node in tip_ids]
    return SeqSim.TreeSequenceSimulation(node_sequences, tip_alignment, root, tip_ids)
end

function SeqSim.simulate_tree_sequences(tree::TreeSim.Tree, site_model::SeqSim.SiteModel)
    return SeqSim.simulate_tree_sequences(Random.GLOBAL_RNG, tree, site_model)
end

"""
    SeqSim.simulate_alignment(rng, tree::TreeSim.Tree, site_model::SeqSim.SiteModel)
    SeqSim.simulate_alignment(tree::TreeSim.Tree, site_model::SeqSim.SiteModel)

Simulate along `tree` and return only the derived tip alignment.

The returned `Alignment` is ordered like `TreeSim.tips(tree)`. Sequence taxon
and time metadata follow the same label fallback and node-time passthrough
policy as `simulate_tree_sequences`.
"""
function SeqSim.simulate_alignment(
    rng::AbstractRNG,
    tree::TreeSim.Tree,
    site_model::SeqSim.SiteModel,
)::SeqSim.Alignment
    return SeqSim.simulate_tree_sequences(rng, tree, site_model).tip_alignment
end

function SeqSim.simulate_alignment(tree::TreeSim.Tree, site_model::SeqSim.SiteModel)
    return SeqSim.simulate_alignment(Random.GLOBAL_RNG, tree, site_model)
end

function SeqSim.simulate_on_tree(
    rng::AbstractRNG,
    tree::TreeSim.Tree,
    site_model::SeqSim.SiteModel,
)::SeqSim.TreeSequenceSimulation
    return SeqSim.simulate_tree_sequences(rng, tree, site_model)
end

function SeqSim.simulate_on_tree(tree::TreeSim.Tree, site_model::SeqSim.SiteModel)
    return SeqSim.simulate_tree_sequences(Random.GLOBAL_RNG, tree, site_model)
end

end
