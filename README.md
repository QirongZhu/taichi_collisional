# taichi_collisional
An N-body code for star cluster simulations. The code uses an adaptive (octree) Fast Multipole Method to compute gravity, therefore $\mathcal{O}(N)$ instead of pairwise summations ($N^2$). The multipole-to-local (M2L) translations uses a $\mathcal{O}(p^3)$ method by rotating the multipoles in the direction of the seperation vector. The particle-to-particle kernel uses a vectorized version. For more details, see (Mukherjee 2021 ApJ, arxiv:2012.02207).

-Much of the code is based on exafmm (https://github.com/exafmm/exafmm). Solid spherical harmonics is adopted over the original formulation (See Dehnen 2014). 

-The master folder contains the code used in the ApJ paper.

-The speed folder includes a fourth-order integrator (arxiv:2011.14984), which improves the energy conservation. The gradient force is approximated using the extraplation method of Omelyan 2006. See also Farr2007. 



