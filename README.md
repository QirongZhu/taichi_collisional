# taichi_collisional
This is a version of parallel FMM aiming for accurate N-body integration.
For a distributed octree, we will use

(0) A Hilbert code hashed octree for a global ID system

(1) A top-tree for work-load partition as well as initial force magnitude estimate;

(2) Given the estimated force, walk the top-tree to find the remote nodes need to be exported;

(3) Communicate the cells and particles to build a locally essential tree.

(4) The rest follows similary as a sequential FMM. 
