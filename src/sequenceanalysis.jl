
"""
    tip_sequences(tree::RootedTree)::Dict{String, Vector{Char}}

Extract sequences from the tips (leaves) of a phylogenetic tree and return them as a dictionary.

# Arguments
- `tree`: A RootedTree object from which tip sequences are to be extracted.

# Returns
- A dictionary mapping tip names to their sequences.

# Example
```julia
tree = ... # some RootedTree with sequences in the node data
tip_seq = tip_sequences(tree)
println(tip_seq)
```
"""
function tip_sequences(tree::AbstractTree)::Dict{String, <: BioSequence} 
    seq_type = typeof(getroot(tree).data["sequence"])   # TODO: Replace this with more robust option case tree doesn't have a root
    seqs = Dict{String, seq_type}()
    for leaf in getleaves(tree)
        if haskey(leaf.data, "sequence")
            seqs[leaf.name] = leaf.data["sequence"]
        else
            @warn "Leaf $(leaf.name) does not contain sequence data."
            seqs[leaf.name] = seq_type()
        end
    end
    return seqs
end