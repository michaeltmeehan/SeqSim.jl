
using Phylo  # Make sure to include this if it's not already in the file

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
function tip_sequences(tree::RootedTree)::Dict{String, Vector{Char}}
    seqs = Dict{String, Vector{Char}}()
    for leaf in getleaves(tree)
        if haskey(leaf.data, "sequence")
            seqs[leaf.name] = leaf.data["sequence"]
        else
            @warn "Leaf $(leaf.name) does not contain sequence data."
            seqs[leaf.name] = Vector{Char}()
        end
    end
    return seqs
end