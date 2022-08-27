#include "tree.h"

int main()
{
    Tree kdtree;
    kdtree.setType(1);
    kdtree.setBodies(1000000);
    
    auto start = std::chrono::steady_clock::now();
    kdtree.buildTree();
    auto stop  =std::chrono::steady_clock::now();
    
    std::chrono::duration<double> e_seconds = stop - start;
    std::cout << "Tree construction took: " << e_seconds.count() << std::endl;
    
    
  return 0;
}
