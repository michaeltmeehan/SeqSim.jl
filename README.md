# SeqSim.jl

`SeqSim.jl` is a Julia package developed for simulating DNA sequences for phylogenetic analysis. It supports generating sequences based on various substitution models, including Jukes-Cantor, Kimura, and more sophisticated models like GTR. The package also integrates sequence evolution along phylogenetic trees, taking into account different rates of evolution and other genomic factors.

## Features
- Supports multiple substitution models (JC, K2P, HKY, GTR, etc.).
- Simulates sequence evolution on phylogenetic trees.
- Provides tools for both stochastic and deterministic simulation approaches.
- Customizable for different evolutionary scenarios.
- Variable site models with gamma distribution for rate heterogeneity.
- Handles invariant sites in sequence simulations.

## Installation

Install `SeqSim.jl` using the Julia package manager. From the Julia REPL, type the following command:

    julia> using Pkg
    julia> Pkg.add("SeqSim")

Alternatively, you can install the latest version directly from GitHub:

    julia> Pkg.add(url="https://github.com/AITHM/SeqSim.jl.git")

## Usage

### Basic Example

Here is a basic example of how to simulate sequences using `SeqSim.jl`:

```
    using Phylo
    using SeqSim

    # Define a substitution model
    model = JC()  # Jukes-Cantor model

    # Create a random phylogenetic tree
    tree = rand(Nonultrametric(10))  # 10 taxa

    # Simulate sequences
    simulate_sequences!(tree, 100, model)  # 100 base pairs long sequences

    # Print the resulting sequences
    sequences = tip_sequences(tree)
    println(sequences)
```

### Advanced Example with Variable Sites

This example shows how to simulate sequences with a model that includes variable sites and rate heterogeneity:
```
    using Phylo
    using SeqSim
    using BioSequences

    # Define a substitution model
    substitution_model = HKY(κ=2.0, π=[0.1, 0.2, 0.3, 0.4])  # Hasegawa-Kishino-Yano model

    # Define a site model with gamma distribution for rate heterogeneity
    site_model = SiteModel(
        mutation_rate = 1.0,
        gamma_category_count = 4,
        gamma_shape = 0.5,
        proportion_invariant = 0.1,
        substitution_model = substitution_model
    )

    # Create a random phylogenetic tree
    tree = rand(Nonultrametric(10))  # 10 taxa

    # Simulate sequences
    simulate_sequences!(tree, 100, site_model)  # 100 base pairs long sequences

    # Print the resulting sequences
    sequences = tip_sequences(tree)
    println(sequences)
```

## Substitution Models

`SeqSim.jl` supports the following substitution models:
- **JC (Jukes-Cantor)**
- **K2P (Kimura 2-Parameter)**
- **HKY (Hasegawa-Kishino-Yano)**
- **GTR (General Time Reversible)**

Each model can be customized with specific evolutionary parameters, such as transition/transversion rates, nucleotide frequencies, and more.

## Variable Site Models

The package supports models with variable sites using gamma distributions for rate heterogeneity and allows the inclusion of invariant sites. This can be configured through the `SiteModel` struct, enabling more realistic simulations of sequence evolution.

## Testing

To run tests after installation, navigate to the package directory and run:

```
    using Pkg
    Pkg.test("SeqSim")
```


## Contributing

Contributions to `SeqSim.jl` are welcome. To contribute:
1. Fork the repository.
2. Create a new branch for your feature.
3. Add your feature or enhancement.
4. Write or update tests as necessary.
5. Submit a pull request.

Please make sure to update tests as appropriate.
