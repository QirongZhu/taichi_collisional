#include "hilbert.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>

int main() 
{
    using namespace hilbert;

    MPI_Init(NULL, NULL);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    std::cout << my_rank << " " << p << std::endl;

    std::vector<double> coords {0.23, 0.5, 0.2};

    int level = 20;
    auto ix = get3DIndex(coords, level);

    std::cout << ix[0] << ix[1] << ix[2] << std::endl;

    MPI_Finalize();

    return 0;
}