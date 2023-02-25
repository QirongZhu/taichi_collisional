#include "hilbert.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <random>
#include <bitset>
#include "allvars.h"
#include "tree.h"

double getRandom() {
    static std::default_random_engine e;
    static std::uniform_real_distribution<> dist(0, 1);
    return dist(e);
}

int main(int argc, char* argv[]) 
{
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(my_processor.my_rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(my_processor.total_ranks));
    hilbert::computeOffsets();
  }

if (my_processor.my_rank == 0)
{
  std::cout << "Welcome to Taichi! This is a parallel version. " << std::endl;
}

  Particles particles(100000);

  uint64_t index = 0;

  for (auto&p : particles)
  {
    p.X[0] = getRandom();
    p.X[1] = getRandom();
    p.X[2] = getRandom();
    p.id= (index++);
    std::array<double, 3> coords{p.X[0], p.X[1], p.X[2]};
    p.hp_key= hilbert::getKey(coords);
    p.done=false;
    if (my_processor.my_rank == -1)
    {
        std::cout << p.hp_key <<" ==> "<< std::bitset<64>(p.hp_key)<< std::endl;
    }
  }  

  topTree top_tree;

  top_tree.buildTopTree(particles);
 
  MPI_Finalize();

  return 0;
}
