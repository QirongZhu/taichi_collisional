#include "tree.h"

using namespace FMM;

int main()
{
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(-10, 10);

  Bodies inputBodies;
  for (int i=0; i<1000000; i++)
  {
    inputBodies.emplace_back(Body(dist(e2), dist(e2), dist(e2), 1.0, i));
  }

  // static method called to initialize some coeff arrays
  Tree::initKernel();

  Tree kdtree;
  kdtree.setType(1);
  kdtree.setBodies(inputBodies);

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
