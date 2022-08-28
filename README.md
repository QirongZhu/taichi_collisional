## Purpose
This is a side-project to investigate which tree format is the best for Fast
Multipole Method (FMM) for high accuracy N-body integration. Historically,
octree is the default underlying data structure for adaptive FMM. The reasons
behind octrees are at least: (1) quick construction time given a list of particles (2) cells/nodes have nice geometries, ie, cubes, s.t. the errors behave as expected from mathmatical analysis.

The downside of octree for FMM is due to the fact that leaf cells have very un-even particle counts. This is not friendly for modern CPUs, which are often memory/bandwidth limited in practice. This leads to load imbalance, and not fully utilize the CPU power. Hence, we set out to explore other options of trees with the goal to find the best data structure for fast N-body solver.

We will consider the following well-known *binary* trees:

### K-d tree

https://en.wikipedia.org/wiki/K-d_tree

### Variant of K-d tree, bisector at center of mass
https://academic.oup.com/mnras/article/418/2/770/1068566

### Press tree (mutual nearest neighour tree)

https://ui.adsabs.harvard.edu/abs/1989ApJS...71..871J/abstract

### Binary radix tree

https://developer.nvidia.com/blog/parallelforall/wp-content/uploads/2012/11/karras2012hpg_paper.pdf

https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/

### Usage
e.g.,

1)Edit CMakeList.txt with "BOOST_DIR" and expansion order "EXPANSION"

2)mkdir build

3)cd build

4)CC=gcc-11 CXX=g++-11 cmake ../ -DCMAKE_BUILD_TYPE=Release

5)make

6)./tree