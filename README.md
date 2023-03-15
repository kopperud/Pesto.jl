## Pesto.jl: Phylogenetic Estimation of Shifts in the Tempo of Origination

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kopperud.github.io/Pesto.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kopperud.github.io/Pesto.jl/dev)

This Julia module also has an R-package wrapper with the same name (see [Pesto](https://github.com/kopperud/Pesto)).

The program can fit state-dependent speciation and extinction (SSE) models on reconstructed phylogenetic trees. We use these models to make inferences about when, and on which branches, there were shifts in the tempo of the process of diversification. The method is equivalent to the one presented by Höhna et al (2019), however it is much faster. The runtime of Pesto scales linearly with the number of tips in the phylogeny, and can be run on huge phylogenies (> 30k taxa) without much trouble on a standard laptop computer.

See the website at [https://kopperud.github.io/Pesto.jl/dev](https://kopperud.github.io/Pesto.jl/dev) for documentation.

## References

* Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.
* Höhna, S., Freyman, W. A., Nolen, Z., Huelsenbeck, J. P., May, M. R., & Moore, B. R. (2019). A Bayesian approach for estimating branch-specific speciation and extinction rates. BioRxiv, 555805.

