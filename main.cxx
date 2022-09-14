#ifdef USE_RADXTREE
#include "radixtree.h"
#else
#include "tree.h"
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <cstring>
#include <cstdlib>


using namespace FMM;

int main(int argc, char **argv)
{
    //to read from a file
    /*static char input_fname[200];
    sprintf(input_fname, argv[1]);


    //read file
    size_t numBodies = 0;
    std::vector < double >array;
    {
        size_t number_of_lines = 0;
        std::string line;
        std::ifstream file(input_fname);
        while(std::getline(file, line))
            ++number_of_lines;
        numBodies = number_of_lines;
        file.close();

        std::ifstream fin(input_fname);
        array.assign(std::istream_iterator < double >(fin), std::istream_iterator < double >());
        fin.close();
    }
*/
    std::cout<<numBodies<<std::endl;
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(0, 1);
    
  Bodies inputBodies;
  for (int i = 0; i < numBodies; i++)
  {
        inputBodies.emplace_back(Body(dist(e2), dist(e2), dist(e2), 1.0, i));
      //inputBodies.emplace_back(Body(array[i * 7+1],array[i * 7+2],array[i * 7+3],array[i * 7+0],i));
  }

    // static method called to initialize some coeff arrays
    //Tree::initKernel();

    Tree kdtree;
    Tree::initKernel();
    kdtree.setType(2);
    kdtree.setBodies(inputBodies);

    auto start = std::chrono::steady_clock::now();
    kdtree.buildTree();
    auto stop = std::chrono::steady_clock::now();
    
    std::chrono::duration<double> e_seconds = stop - start;
    std::chrono::duration<double> e_total = e_seconds;
    std::cout << "Tree construction took: " << e_seconds.count() << std::endl;

    start = std::chrono::steady_clock::now();
    kdtree.upwardPass();
    stop = std::chrono::steady_clock::now();
    e_seconds = stop - start;
    e_total += e_seconds;
    std::cout << "Tree upward traversal took: " << e_seconds.count() << std::endl;

    // kdtree.printTree();

    start = std::chrono::steady_clock::now();
    kdtree.horizontalPass();
    stop = std::chrono::steady_clock::now();
    e_seconds = stop - start;
    e_total += e_seconds;
    std::cout << "Tree horizontal traversal took: " << e_seconds.count() << std::endl;
    
    start = std::chrono::steady_clock::now();
    kdtree.downwardPass();
    stop = std::chrono::steady_clock::now();
    e_seconds = stop - start;
    e_total += e_seconds;
    std::cout << "Tree downward traversal took: " << e_seconds.count() << std::endl;


    std::cout<<"Total time:" << e_total.count() <<std::endl;
    //kdtree.dumpForces(numBodies); 

    // kdtree.printTree();
  
  return 0;
}
