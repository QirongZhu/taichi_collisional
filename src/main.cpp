#include "hilbert.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include "allvars.h"

Domain domain;

struct myProcessor
{
  int my_rank;
  int total_ranks;
} my_processor;


int main() 
{
  using namespace hilbert;

  {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &(my_processor.my_rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(my_processor.total_ranks));
    std::cout << my_processor.my_rank << " out of " << my_processor.total_ranks << std::endl;
  }

  MPI_Finalize();

  return 0;
}
