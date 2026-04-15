module SeqSim

using Crayons
using Distributions
using LinearAlgebra
using Random
using StaticArrays


include("sequence.jl")
include("processing.jl")
include("inspection.jl")
include("viewer.jl")
include("substitution.jl")
include("site.jl")
include("sequence_simulator.jl")
include("SequencePropagator.jl")
include("export.jl")


export JC, F81, K2P, HKY, GTR, SubstitutionModel
export SiteModel
export rand_seq, SequencePropagator
export Sequence
export filter_sequences, filter_sites
export variable_sites, invariant_sites, get_snps, snp_alignment, invariant_filtered_alignment
export site_states, site_snp_counts, sequence_snp_counts
export site_state_counts, observed_state_counts
export minor_allele_counts, minor_allele_frequencies, nonreference_frequencies
export singleton_sites, parsimony_informative_sites
export pairwise_differences, pairwise_difference_matrix, consensus_sequence
export AlignmentSummary
export site_window, centered_site_window
export SiteInspectionRow, AlignmentInspection, inspect_sites, inspect_window, inspect_centered_window
export selected_site_strings
export AlignmentViewer, alignment_viewer, window_viewer, centered_window_viewer, variable_site_viewer
export viewer_preview
export write_fasta, write_nexus, write_phylip, write_alignment

end # module SeqSim
