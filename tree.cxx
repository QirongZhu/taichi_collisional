#include "tree.h"

using namespace FMM;

int main()
{

  // static method called to initialize some coeff arrays
  Tree::initKernel();

  Tree kdtree;
  kdtree.setType(1);
  kdtree.setBodies(1000000);

  auto start = std::chrono::steady_clock::now();
  kdtree.buildTree();
  auto stop = std::chrono::steady_clock::now();

  std::chrono::duration<double> e_seconds = stop - start;
  std::cout << "Tree construction took: " << e_seconds.count() << std::endl;

  start = std::chrono::steady_clock::now();
  kdtree.upwardPass();
  stop = std::chrono::steady_clock::now();
  e_seconds = stop - start;
  std::cout << "Tree upward traversal took: " << e_seconds.count() << std::endl;

  kdtree.printTree();

  return 0;
}
