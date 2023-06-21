# Home

## Background 
Pesto is an acronym for Phylogenetic Estimation of Shifts in the Tempo of Origination. 
Broadly speaking, it is a method for detecting shifts in the process of diversification that led to the biodiversity present today.
In order to study diversification, we are interested in two events: speciation and extinction events. 
Since these events are difficult to observe, we use a simple model, the birth-death process, to model what happened in the past. 
The tempo at which species speciate and go extinct are controlled by two parameters (Nee et al. 1994):
* The speciation rate ($\lambda$)
* The extinction rate ($\mu$)
By simulating under the birth-death process, and pruning the extinct lineages from the tree, one gets a **reconstructed phylogenetic tree** as a result (i.e. it is ultrametric, all tips end at the same time point).

The question we are interested in, is whether the process of diversification changed throughout the phylogenetic tree. 
In other words, was there a shift or not, and if so, how large was the shift?
To do so, we are employing a variant of the state-dependent birth-death model (first presented by Maddison et al. 2007).
This is also called the birth-death-shift model, or the lineage-specific birth-death model (Höhna et al. 2019).
The birth-death-shift model has an additional parameter:
* The shift rate ($\eta$). 
When a shift event occurs, the speciation and extinction rate shifts from the previous state (say $\lambda_1,\mu_1$) to a new state with different rates ($\lambda_2,\mu_2$).

## References

* Nee, S., May, R. M., & Harvey, P. H. (1994). The reconstructed evolutionary process. Philosophical Transactions of the Royal Society of London. Series B: Biological Sciences, 344(1309), 305-311.
* Maddison, W. P., Midford, P. E., & Otto, S. P. (2007). Estimating a binary character's effect on speciation and extinction. Systematic biology, 56(5), 701-710.
* Höhna, S., Freyman, W. A., Nolen, Z., Huelsenbeck, J. P., May, M. R., & Moore, B. R. (2019). A Bayesian approach for estimating branch-specific speciation and extinction rates. BioRxiv, 555805.