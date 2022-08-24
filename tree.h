#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

struct Particle
{
  double X[3];
  double M;
};

typedef std::vector<Particle> Particles;

struct Cell
{
  Cell* CHILD;
  int  NCHILD;
  Particle* BODY;
  int  NBODY;
};

typedef std::vector<Cell> Cells;
