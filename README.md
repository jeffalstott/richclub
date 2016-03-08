# richclub
Python package for rich club analysis and detection of graphs, using igraph

This is a smattering of code that was used in the process of creating the paper:

Jeff Alstott, Pietro Panzarasa, Mikail Rubinov, Edward T. Bullmore & Petra E. Vértes. (2014). A Unifying Framework for Measuring Weighted Rich Clubs. [Scientific Reports 4: 7258](http://www.nature.com/articles/srep07258)

(Academics, if you use this code please cite the paper.)

This code is really overgrown and is mainly an artifact of my learning process as we figured out how to best measure rich clubs. Now that the learning process is done, all that remains is very simple:

1. Measure the amount of weight on the links between nodes of a particular strength rank or higher. The sum of these weights is called `phi`.
2. Consider what the randomized controls should be, and make a lot of randomized controls. Calculate `phi` for each control, and average it.
3. Divide the `phi` of the first network in step #1 to the average `phi` calculated from the controls in step #2. This is the normalized rich club coefficient.

Code for calculating rich clubs
===

Simple case
---
If you have a network represented a `numpy` array `X` and you want to calculate the `phi` at every possible strength threshold, just use `richclub.fast_RC(X)`.

Complex case
---
If you have a network represented as an `igraph` network `X`, and you want to calculate the rich club structure, you can instead create an object for analysis with `this_analysis = richclub.RC(X)`. You can then use `this_analysis.phis()` to get the sum of the weights of the links between nodes of a strength rank or higher, for all possible ranks. But `richclub.RC` also allows you customize. Instead of using strengths, you can give each node a precalculated score to threshold on instead: `this_analysis = richclub.RC(X, scores=range(10))`. You can also specificy which threshold ranks to use (instead of just every unique value of score) by specifying them: `this_analysis = richclub.RC(X, scores=range(10), ranks=range(5))`

Randomized controls
---
The hard part of calculating rich clubs is determining what is an appropriate randomized control, and then calculating it. The [paper](http://www.nature.com/articles/srep07258) talks about this at length. If you have a network `X`, you can use `richclub.preserve_strength(X)` to generate randomized versions of `X` using a variety of different options. It is rather heavy-duty. Options include whether the controls should:
- randomize topology (`randomize_topology=True`, randomization method set for `igraph` with `randomize_method='vl'`)
- permute the strengths of the nodes among nodes with the same degrees (`permute_strength=True`)
- exactly preserve the out-strengths or the in-strengths of the original network (`preserve_mode='out'` vs. `preserve_mode='in'`) or instead approximately preserve both (`preserve_mode='estimate_both'`, method from [Serrano et al. 2006, equation 10](http://arxiv.org/abs/cond-mat/0609029)).


What a basic workflow could look like
-----
```
from richclub import preserve_strength, RC
from numpy import mean
empirical_phis = RC(X, scores=X.degree()).phis()
average_randomized_phis = mean([RC(preserve_strength(X), scores=X.degree()).phis() for i in range(1000)], axis=0)
normalized_phis = empirical_phis/average_randomized_phis
```

This would calculate the weighted rich clubs of a network, thresholded by degree, normalized by the randomized case where the topology and strength sequence are preserved but weights are allowed to vary. 


Dependencies
===
- numpy
- igraph


