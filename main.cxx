//#include "tree.h"
#include "radixtree.h"

using namespace FMM;

int main()
{
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(0, 1);

  Bodies inputBodies;
  for (int i = 0; i < 1000; i++)
  {
    inputBodies.emplace_back(Body(dist(e2), dist(e2), dist(e2), 1.0, i));
  }

  Tree bintree;
  bintree.setBodies(inputBodies);
  //  bintree.getBoundBox();
  auto start = std::chrono::steady_clock::now();

  bintree.buildRadixTree();
  std::cout << "buildRadixTree done\n";

  auto stop = std::chrono::steady_clock::now();
  std::chrono::duration<double> e_seconds = stop - start;
  std::cout << "Tree construction took: " << e_seconds.count() << std::endl;

  // bintree.printTree();
  start = std::chrono::steady_clock::now();

  bintree.sumUpward();
  std::cout << "upwardPass done\n";

  bintree.flagNode();
  std::cout << "flag big nodes \n\n ";
  stop = std::chrono::steady_clock::now();

  e_seconds = stop - start;
  std::cout << "Tree postprocessing took: " << e_seconds.count() << std::endl;

  bintree.convertCells();

  /*
  #ifdef USE_OCTREE
    std::cout << "Octree used for FMM \n";
  #else
    std::cout << "Binary used for FMM \n";
  #endif

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

    // kdtree.printTree();

    start = std::chrono::steady_clock::now();
    kdtree.horizontalPass();
    stop = std::chrono::steady_clock::now();
    e_seconds = stop - start;
    std::cout << "Tree horizontal traversal took: " << e_seconds.count() << std::endl;

    // kdtree.printTree();
  */
  return 0;
}
