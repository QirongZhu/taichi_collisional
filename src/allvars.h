#pragma once
#include <vector>
#include <numeric>

/*
  definitions of structs
*/

const double MAXDOUBLE = std::numeric_limits<double>::max();


struct myProcessor
{
  int my_rank;
  int total_ranks;
};

extern myProcessor my_processor;

struct Domain {
  double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
  double X0[3];
  double R0;
};


struct Particle {
    double X[3];
    uint64_t id;
    uint64_t hp_key;
    uint64_t top_node;
    bool done=false;
};

typedef std::vector<Particle> Particles;

extern Domain domain;
