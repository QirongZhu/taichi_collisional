#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <cstring>

#include "exafmm.h"
#include "integrator.h"
#include "kernel_fast_lean.h"
#include "timer.h"
#include "traverse_eager.h"
#include "load_gadget.h"
#include "io.h"

using namespace std;
using namespace exafmm;

void init_code(void);

double interval;

struct diagnostics global_diag;
struct diagnostics *diag;
struct sys zerosys = { 0, NULL, NULL };

int main(int argc, char **argv)
{
  P = EXPANSION;

  long double cmpos[3], cmvel[3];

  initKernel();
  init_code();

  printf("\n");
  printf("This is Taichi_cosmo v0.5. Taichi_cosmo is a FMM based N-body code with individual timesteps.\n");
  printf("Author: Qirong Zhu.\n");
  printf("    ;:::      \n");
  printf("   #+    :   \n");
  printf("  @#       `  \n");
  printf(" ###  ##    : \n");
  printf(" ###  .,       \n");
  printf(":###.        :\n");
  printf("#####,       :\n");
  printf("########+   :\n");
  printf(":########@   :\n");
  printf(" ######+###    \n");
  printf(" #####   ## : \n");
  printf("  @####:## `  \n");
  printf("   ######.:   \n");
  printf("     :##,     \n");
  printf("\n");
  fflush(stdout);

  if(argc != 5)
    {
      printf("\nwrong number of arguments\n");
      printf
	("call with:\n\n ./Taichi 0 <initial_condition_name> <t_end> <output_interval> for Gadget  snapshot or \n");
      printf("               ./Taichi 1 <snapshot_name> <t_end> <output_interval> for restart snapshot \n\n");
      printf("number of arguments %d\n", argc);
      fflush(stdout);

      exit(0);
    }

  sprintf(input_fname, argv[2]);

  double t_end = std::stof(argv[3]);
  interval = std::stof(argv[4]);
  double dt = interval;

  if(std::stoi(argv[1]) == 0)
    {
      snapnum = 0;
      t_now = 0;

      bool gadget_file = false;
      std::vector < double >array;

      if(!gadget_file)
	{
	  start("reading bodies");
	  size_t number_of_lines = 0;
	  std::string line;
	  std::ifstream file(input_fname);
	  while(std::getline(file, line))
	    ++number_of_lines;
	  numBodies = number_of_lines;
	  file.close();

	  std::ifstream fin(input_fname);
	  array.resize(numBodies * 7);
	  array.assign(std::istream_iterator < double >(fin), std::istream_iterator < double >());
	  fin.close();
	  stop("reading bodies");
	}
      else
	{
	  int files = 1;
	  load_snapshot(input_fname, files);
	  printf("Reading IC done with %d particles \n", NumPart);
	  fflush(stdout);
	  numBodies = NumPart;
	}

      mainsys.n = numBodies;
      mainsys.part = (struct particle *) malloc(numBodies * sizeof(struct particle));
      mainsys.last = &mainsys.part[numBodies - 1];

      if(!gadget_file)
	{
	  for(size_t b = 0; b < numBodies; b++)
	    {
	      mainsys.part[b].id = b;
	      mainsys.part[b].mass = array[b * 7 + 0];
	      mainsys.part[b].pos[0] = array[b * 7 + 1];
	      mainsys.part[b].pos[1] = array[b * 7 + 2];
	      mainsys.part[b].pos[2] = array[b * 7 + 3];
	      mainsys.part[b].vel[0] = array[b * 7 + 4];
	      mainsys.part[b].vel[1] = array[b * 7 + 5];
	      mainsys.part[b].vel[2] = array[b * 7 + 6];
	    }
	}
      else
	{
	  for(size_t b = 0; b < numBodies; b++)
	    {
	      mainsys.part[b].pos[0] = PPP[b + 1].Pos[0];
	      mainsys.part[b].pos[1] = PPP[b + 1].Pos[1];
	      mainsys.part[b].pos[2] = PPP[b + 1].Pos[2];
	      mainsys.part[b].vel[0] = PPP[b + 1].Vel[0];
	      mainsys.part[b].vel[1] = PPP[b + 1].Vel[1];
	      mainsys.part[b].vel[2] = PPP[b + 1].Vel[2];
	      mainsys.part[b].mass = PPP[b + 1].Mass;
	      mainsys.part[b].id = PPP[b + 1].Id;
	    }
	}
    }
  else
    {
      read_snapshot(snapnum);
    }

  for(unsigned int b = 0; b < mainsys.n; b++)
    {
      mainsys.part[b].pos_e[0] = 0;
      mainsys.part[b].pos_e[1] = 0;
      mainsys.part[b].pos_e[2] = 0;
    }

  start("Intergration");

  if(t_end > 0)
    {
      while(t_end > t_now)
	{

#ifdef DEBUG
	  start("Dummy Poisson test");
	  kick_naive(0, mainsys, zerosys, zerosys, 0, false);
	  stop("Dummy Poisson test");
	  fflush(stdout);
#endif	  

	  write_snapshot(snapnum, mainsys);
	  snapnum++;

	  real_t kinetic = system_kinetic_energy(mainsys);
	  real_t pot = system_potential_energy(mainsys);
	  system_center_of_mass(mainsys, cmpos, cmvel);
	  printf("pot=%18.12f kin=%18.12f tot=%18.12f \n", pot, kinetic, -pot + kinetic);
	  fflush(stdout);

	  kick_naive(0, mainsys, zerosys, zerosys, 0, true);

	  if(mainsys.n > 0)
	    do_evolve(mainsys, dt);

	  t_now += dt;
	}
    }

  stop("Intergration");

  return 0;
}

void init_code(void)
{
  diag = &global_diag;
  diag->simtime = 0.;
}
