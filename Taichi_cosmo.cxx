#include "build_tree.h"
#include "new_kernel_cart.h"
#include "timer.h"
#include "load_gadget.h"
#include "evolve.h"
#include "hdf5.h"
#include "traverse_eager.h"

#include <iostream>
#include <sys/stat.h>

using namespace std;
using namespace exafmm;

struct diagnostics global_diag;
struct diagnostics *diag;
struct sys zerosys = { 0, -HUGE, HUGE, NULL, NULL };

double t_now;
double interval;

static int snapnum, snapnum0;
static char input_fname[200];
static int number_files_to_read;
static int number_files_to_write;

void init_code(void);
void init_param(void);
void system_center_of_mass(struct sys s, long double *cmpos, long double *cmvel);
real_t mean_velocity(struct sys s);
void write_snapshot(int snapshot_num, struct sys s, struct io_header header);

int main(int argc, char **argv)
{
  printf("\n");
  printf("This is Taichi_cosmo v0.5. Taichi_cosmo is a FMM based N-body code with individual timesteps.\n");
  printf("Author: Qirong Zhu.\n");
  printf("    ;:::      \n");
  printf("   #+     :   \n");
  printf("  @#       `  \n");
  printf(" ###  ##    : \n");
  printf(" ###  .,       \n");
  printf(":###.        :\n");
  printf("#####,       :\n");
  printf("########+    :\n");
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

  unsigned int numBodies;
  struct sys mainsys;
  long double cmpos[3], cmvel[3];

  srand48(124);			// Set seed for random number generator

  init_param();

  snapnum = 0;
  snapnum0 = 0;

  if(std::stoi(argv[1]) == 0)
    {
      load_snapshot(input_fname, number_files_to_read);

      printf("Reading IC done\n");
      printf("total numbers of particles =%d at t=%g \n", NumPart, Time);
      printf("box=%g omega_0=%g omega_l=%g \n\n", header.BoxSize, header.Omega0, header.OmegaLambda);
      //      printf("masstable=%g %g %g \n", header.mass[0], header.mass[1], header.mass[2]);
      //      printf("pos=%g %g %g mass=%g \n", PPP[0].Pos[0],  PPP[0].Pos[1],  PPP[0].Pos[2], PPP[0].Mass);
      //      printf("pos=%g %g %g mass=%g \n", PPP[NumPart-1].Pos[0],  PPP[NumPart-1].Pos[1],  PPP[NumPart-1].Pos[2], PPP[NumPart-1].Mass);
      fflush(stdout);

      cycle = (real_t) header.BoxSize;

      numBodies = NumPart;

      totnumBodies = NumPart;

      initKernel();

      init_code();

      printf("numBodies=%d NumPart=%d\n", numBodies, NumPart);
      printf("particle_data=%ld %d\n", sizeof(particle_data), numBodies);
      printf("particle=%ld %d\n", sizeof(particle), numBodies);
      fflush(stdout);

      printf("particle_data%ld %d\n", sizeof(particle_data), numBodies);
      printf("particle=%ld %d\n", sizeof(particle), numBodies);
      fflush(stdout);

      start("Setup Main system");

      mainsys.n = numBodies;

      mainsys.part = (struct particle *) 
	malloc(numBodies * sizeof(struct particle));

      mainsys.last = &mainsys.part[numBodies - 1];

      real_t vel_fac = sqrt(Time) * Time;
#pragma omp parallel for
      for(size_t b = 0; b < numBodies; b++)
	{
	  mainsys.part[b].pos[0] = (real_t) PPP[b].Pos[0];
	  mainsys.part[b].pos[1] = (real_t) PPP[b].Pos[1];
	  mainsys.part[b].pos[2] = (real_t) PPP[b].Pos[2];

	  mainsys.part[b].vel[0] = (real_t) PPP[b].Vel[0] * vel_fac;
	  mainsys.part[b].vel[1] = (real_t) PPP[b].Vel[1] * vel_fac;
	  mainsys.part[b].vel[2] = (real_t) PPP[b].Vel[2] * vel_fac;

	  mainsys.part[b].id     = PPP[b].ID;

	  //mainsys.part[b].mass   =(real_t)PPP[b].Mass;     
	  //mainsys.part[b].type   =PPP[b].Type; 
	}

      stop("Setup Main system");
      fflush(stdout);

      free_gadget_memory();
    }
  else
    {
#include "read_snapshot.h"
    }

  init_drift_table();

  bodies.reserve(numBodies);
  bodies.resize(numBodies);
  cells.reserve(numBodies / (ncrit / 5));

  masspart = header.mass[1];

  X0[0] = cycle * 0.5;
  X0[1] = cycle * 0.5;
  X0[2] = cycle * 0.5;
  R0 = cycle * 0.5 * 1.0001;

  start("copy bodies");

  for(size_t i = 0; i < numBodies; i++)
    {
      for(int d = 0; d < 3; d++)
	bodies[i].X[d] = mainsys.part[i].pos[d];

      bodies[i].issink = 1;
      bodies[i].index = i;
      //bodies[i].q       = mainsys.part[i].mass;
      //bodies[i].type    = mainsys.part[i].type;
      bodies[i].issource = 1;
    }

  printf("cycle=%g center=%g %g %g numBodies=%d\n", cycle, X0[0], X0[1], X0[2], numBodies);

  stop("copy bodies");
  fflush(stdout);

  get_force_and_potential(bodies, 1, 0);
  fflush(stdout);

  real_t acc_mag, min_dt, max_dt;

  max_dt = -HUGE;
  min_dt = HUGE;

  real_t a2inv = 1.0 / (Time * Time);
  real_t a3inv = a2inv / Time;
  real_t hubble_a = diag->Omega0 * a3inv + 
    diag->OmegaLambda + (1 - diag->Omega0 - diag->OmegaLambda) * a2inv;

  hubble_a = Hubble * sqrt(hubble_a);

#pragma omp parallel for reduction(min:min_dt) reduction(max:max_dt)
  for(size_t b = 0; b < bodies.size(); b++)
    {
      unsigned int i = bodies[b].index;

      acc_mag = sqrt(bodies[b].d.F[0] * bodies[b].d.F[0] +
		     bodies[b].d.F[1] * bodies[b].d.F[1] + 
		     bodies[b].d.F[2] * bodies[b].d.F[2]);

      mainsys.part[i].acc_part[0] = 0;
      mainsys.part[i].acc_part[1] = 0;
      mainsys.part[i].acc_part[2] = 0;

      mainsys.part[i].acc_mag = acc_mag;

      mainsys.part[i].timestep = sqrt(eta * Time * softtype[1] / (acc_mag * a2inv * G));

      mainsys.part[i].timestep *= hubble_a;

      if(mainsys.part[i].timestep > MaxSizeTimestep)
	mainsys.part[i].timestep = MaxSizeTimestep;

      if(mainsys.part[i].timestep < MinSizeTimestep)
	mainsys.part[i].timestep = MinSizeTimestep;

      mainsys.part[i].timestep_all = mainsys.part[i].timestep;
      
      if(mainsys.part[i].timestep < min_dt)
	min_dt = mainsys.part[i].timestep;

      if(mainsys.part[i].timestep > max_dt)
	max_dt = mainsys.part[i].timestep;
    }

  mainsys.min_dt = min_dt;

  mainsys.max_dt = max_dt;

  t_now = diag->simtime;

  double t_fric = 0;
  do{  
    t_fric -= interval;
  }while(t_fric > diag->begintime);

  double dt = fmax(interval, max_dt);

  printf("min dt is %g max dt is %g dt for snap0 is %g \n", 
	 min_dt, max_dt,
	 fmin(dt, t_fric+interval-t_now));
  fflush(stdout);

  if(snapnum==0){
    dt = fmin(dt, t_fric+interval-t_now);        
    do_evolve(mainsys, dt);
    t_now += dt;
  }

  t_end += dt / 1048576;

  printf("t_end=%g t_now=%g snapnum=%d \n", t_end, t_now, snapnum);
  fflush(stdout);

  start("Start intergration");

  while(SIGN(dt) * (t_end - t_now) > SIGN(dt) * dt / 1048576)
    {
      simtime = exp(t_now);

      hubble_a = diag->Omega0 / (simtime * simtime * simtime) +
	diag->OmegaLambda + (1 - diag->Omega0 - diag->OmegaLambda) / (simtime * simtime);

      hubble_a = Hubble * sqrt(hubble_a);

      a2inv = 1.0 / (simtime * simtime);
      a3inv = a2inv / simtime;

      system_center_of_mass(mainsys, cmpos, cmvel);

      printf("t=%-6.2f COM %-12.7Lf %-12.7Lf %-12.7Lf ", diag->simtime, cmpos[0], cmpos[1], cmpos[2]);
      printf("VEL %-12.7Lf %-12.7Lf %-12.7Lf \n", cmvel[0], cmvel[1], cmvel[2]);
      fflush(stdout);


      for(size_t i = 0; i < numBodies; i++)
	{
	  for(int d = 0; d < 3; d++)
	    bodies[i].X[d] = mainsys.part[i].pos[d];
	  bodies[i].issink = 1;
	  bodies[i].index = i;
	  //bodies[i].q       = mainsys.part[i].mass;
	  //bodies[i].type    = mainsys.part[i].type;
	  bodies[i].issource = 1;
	}

      get_force_and_potential(bodies, 1, 0);

      max_dt = -HUGE;
      min_dt =  HUGE;

#pragma omp parallel for reduction(min:min_dt) reduction(max:max_dt)
      for(size_t b = 0; b < bodies.size(); b++)
	{
	  unsigned int i = bodies[b].index;

	  acc_mag = sqrt(bodies[b].d.F[0] * bodies[b].d.F[0] +
			 bodies[b].d.F[1] * bodies[b].d.F[1] +
			 bodies[b].d.F[2] * bodies[b].d.F[2]);

	  mainsys.part[i].acc_mag = acc_mag;

	  mainsys.part[i].timestep = sqrt(eta * Time * softtype[1] / (acc_mag * a2inv * G));

	  mainsys.part[i].timestep *= hubble_a;

	  mainsys.part[i].acc_part[0] = 0;
	  mainsys.part[i].acc_part[1] = 0;
	  mainsys.part[i].acc_part[2] = 0;

	  if(mainsys.part[i].timestep > MaxSizeTimestep)
	    mainsys.part[i].timestep = MaxSizeTimestep;

	  if(mainsys.part[i].timestep < MinSizeTimestep)
	    mainsys.part[i].timestep = MinSizeTimestep;

	  mainsys.part[i].timestep_all = mainsys.part[i].timestep;

	  if(mainsys.part[i].timestep < min_dt)
	    min_dt = mainsys.part[i].timestep;

	  if(mainsys.part[i].timestep > max_dt)
	    max_dt = mainsys.part[i].timestep;
	}

      mainsys.min_dt = min_dt;

      mainsys.max_dt = max_dt;

      if(t_now > (snapnum - snapnum0) * interval - TINY + diag->begintime)
	{
	  header.time = exp(t_now);
	  header.redshift = 1.0 / exp(t_now) - 1.0;
	  write_snapshot(snapnum, mainsys, header);
	  snapnum++;
	}

      if(mainsys.n > 0)
	{
	  real_t boxsize = diag->BoxSize;

	  shift_vector[0] = (drand48() - 0.5) * 0.05 * cycle;
	  shift_vector[1] = (drand48() - 0.5) * 0.05 * cycle;
	  shift_vector[2] = (drand48() - 0.5) * 0.05 * cycle;
	  printf("\n shifting vector is: %g %g %g \n", shift_vector[0], shift_vector[1], shift_vector[2]);
	  fflush(stdout);

	  for(size_t b = 0; b < mainsys.n; b++)
	    {
	      mainsys.part[b].pos[0] += shift_vector[0];
	      mainsys.part[b].pos[1] += shift_vector[1];
	      mainsys.part[b].pos[2] += shift_vector[2];

	      if(mainsys.part[b].pos[0] < 0)
		mainsys.part[b].pos[0] += boxsize;
	      if(mainsys.part[b].pos[1] < 0)
		mainsys.part[b].pos[1] += boxsize;
	      if(mainsys.part[b].pos[2] < 0)
		mainsys.part[b].pos[2] += boxsize;
	      if(mainsys.part[b].pos[0] > boxsize)
		mainsys.part[b].pos[0] -= boxsize;
	      if(mainsys.part[b].pos[1] > boxsize)
		mainsys.part[b].pos[1] -= boxsize;
	      if(mainsys.part[b].pos[2] > boxsize)
		mainsys.part[b].pos[2] -= boxsize;
	    }

	  do_evolve(mainsys, dt);

	  for(size_t b = 0; b < mainsys.n; b++)
	    {
	      mainsys.part[b].pos[0] -= shift_vector[0];
	      mainsys.part[b].pos[1] -= shift_vector[1];
	      mainsys.part[b].pos[2] -= shift_vector[2];

	      if(mainsys.part[b].pos[0] < 0)
		mainsys.part[b].pos[0] += boxsize;
	      if(mainsys.part[b].pos[1] < 0)
		mainsys.part[b].pos[1] += boxsize;
	      if(mainsys.part[b].pos[2] < 0)
		mainsys.part[b].pos[2] += boxsize;
	      if(mainsys.part[b].pos[0] > boxsize)
		mainsys.part[b].pos[0] -= boxsize;
	      if(mainsys.part[b].pos[1] > boxsize)
		mainsys.part[b].pos[1] -= boxsize;
	      if(mainsys.part[b].pos[2] > boxsize)
		mainsys.part[b].pos[2] -= boxsize;
	    }
	}

      t_now += dt;
    }

  stop("Start intergration");

  free(mainsys.part);

  return 0;
}


void init_param(void)
{
#include "params.h"
  images = 5;
  Eps_inv = 1.0f / ((float) softtype[1] * 3.0f);
  Eps_inv3 = Eps_inv * Eps_inv * Eps_inv;
}

void init_code(void)
{
  diag = &global_diag;
  diag->simtime = 0.;
  diag->Omega0 = header.Omega0;
  diag->OmegaLambda = header.OmegaLambda;
  diag->BoxSize = header.BoxSize;
  diag->begintime = log(Time);
  diag->simtime = log(Time);
  simtime = Time;
}

void system_center_of_mass(struct sys s, long double *cmpos, long double *cmvel)
{
  long double mass = 0., pos[3] = { 0., 0., 0. }, vel[3] =
  {
  0., 0., 0.};

  for(size_t p = 0; p < s.n; p++)
    {
      for(int i = 0; i < 3; i++)
	{
	  pos[i] += (long double) masspart *s.part[p].pos[i];
	  vel[i] += (long double) masspart *s.part[p].vel[i];
	}
      mass += (long double) masspart;
    }
  for(int i = 0; i < 3; i++)
    {
      cmpos[i] = pos[i] / mass;
      cmvel[i] = vel[i] / mass;
    }
}

real_t mean_velocity(struct sys s)
{
  long double count = 0.0, vel = 0.0;

  for(size_t p = 0; p < s.n; p++)
    {
      vel +=
	s.part[p].vel[0] * s.part[p].vel[0] + s.part[p].vel[1] * s.part[p].vel[1] +
	s.part[p].vel[2] * s.part[p].vel[2];
      count += 1.0;
    }
  return (real_t) (sqrt(vel / count));
}

void write_snapshot(int snapshot_num, struct sys s, struct io_header header)
{
#include "write_snapshot.h"
}
