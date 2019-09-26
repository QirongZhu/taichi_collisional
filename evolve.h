#ifndef evolve_h
#define evolve_h

#include "exafmm.h"
#include "traverse_eager.h"

using namespace exafmm;

#define eta  ((real_t) 0.05)
#define MaxSizeTimestep    ((real_t) 0.025)
#define MinSizeTimestep    ((real_t) 1.e-8)
#define TINY ((real_t) 1e-30)
#define HUGE ((real_t) 1.e30)

#define DRIFT_TABLE_LENGTH 1000
#define Hubble 100

extern struct sys zerosys;

struct __attribute__ ((__aligned__(16))) particle
{
  real_t pos[3];
  unsigned int id;
  real_t vel[3];
  float timestep;
  float acc_part[3];
  float acc_mag;
  float timestep_all;  
#if INDI
  int Ti_begstep;
  int Ti_endstep;
  bool isactive;
#endif
};

struct diagnostics
{
  real_t simtime;
  real_t BoxSize;
  real_t Omega0;
  real_t OmegaLambda;
  real_t begintime;
};				/* diagnostics */

extern struct diagnostics global_diag;
extern struct diagnostics *diag;

#pragma omp threadprivate(diag)

struct sys
{
  size_t n;
  real_t max_dt;
  real_t min_dt;
  struct particle *part;
  struct particle *last;
};

void init_drift_table(void);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double get_gravkick_factor(double time0, double time1);
double get_drift_factor(double time0, double time1);

void kick(int clevel, struct sys s1, struct sys s2, real_t dt,
	  int fast_or_slow, int first_or_second, real_t t1, real_t t2);
void drift(int clevel, struct sys s, real_t etime, real_t dt);

void kdk(int clevel, struct sys s1, struct sys s2, real_t stime, real_t etime, real_t dt);
void dkd(int clevel, struct sys s1, struct sys s2, real_t stime, real_t etime, real_t dt);
static void split(real_t dt, struct sys s, struct sys *slow, struct sys *fast);
void evolve_split_hold(int clevel, struct sys s, real_t stime, real_t etime, real_t dt, int calc_timestep);
void evolve_split_hold_dkd(int clevel, struct sys s, real_t stime,
			   real_t etime, real_t dt, int calc_timestep);
void do_evolve(struct sys s, double dt);

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}
#define ABS(X) (((X) >= 0) ? (X) : -(X))
#define SIGN(X)   ((X>0)-(X<0))

void get_force_and_potential(Bodies & bodies, bool verbose, int high_accuracy)
{
  if(!verbose)
    {
      buildTree(bodies);
      upwardPass(cells);

      if(high_accuracy)
	{
	  for(size_t i = 0; i < cells.size(); i++)
	    {
	      for(int n = 0; n <= P; n++)
		{
		  cells[i].Pn[n] = cells[i].CellMass * std::pow(cells[i].R, n);
		}
	    }
	}
      horizontalPass(cells, cells, high_accuracy);
      downwardPass(cells);
    }
  else
    {
      start("FMM");
      start("build_tree");
      buildTree(bodies);
      stop("build_tree");
      start("P2M and M2M");
      upwardPass(cells);
      stop("P2M and M2M");

      if(high_accuracy)
	{
	  for(size_t i = 0; i < cells.size(); i++)
	    {
	      for(int n = 0; n <= P; n++)
		{
		  cells[i].Pn[n] = cells[i].CellMass * std::pow(cells[i].R, n);
		}
	    }
	}

      start("M2L and P2P");
      horizontalPass(cells, cells, high_accuracy);
      stop("M2L and P2P");
      start("L2L and L2P");
      downwardPass(cells);
      stop("L2L and L2P");
      stop("FMM");
    }

  real_t monopole = 0.0;
  real_t dipole_x = 0.0, dipole_y = 0.0, dipole_z = 0.0;

#pragma omp parallel for reduction(+:monopole,dipole_x,dipole_y,dipole_z)
  for(size_t b = 0; b < bodies.size(); b++)
    {
      monopole += bodies[b].issource;
      dipole_x += bodies[b].X[0] * bodies[b].issource;
      dipole_y += bodies[b].X[1] * bodies[b].issource;
      dipole_z += bodies[b].X[2] * bodies[b].issource;
    }

  real_t coef = 4.0 * M_PI / (3.0 * cycle * cycle * cycle);

  monopole *= coef * masspart;

  dipole_x *= masspart;
  dipole_y *= masspart;
  dipole_z *= masspart;

#pragma omp parallel for
  for(size_t b = 0; b < bodies.size(); b++)
    {      
      bodies[b].d.F[0] -= coef * dipole_x;
      bodies[b].d.F[1] -= coef * dipole_y;
      bodies[b].d.F[2] -= coef * dipole_z;

      bodies[b].d.F[0] += monopole * bodies[b].X[0];
      bodies[b].d.F[1] += monopole * bodies[b].X[1];
      bodies[b].d.F[2] += monopole * bodies[b].X[2];
    }
}

void kdk(int clevel, struct sys s1, struct sys s2, real_t stime, real_t etime, real_t dt)
{
  int first_or_second = 0;
  int fast_or_slow = 1;

  if(s2.n > 0)
    {
      fast_or_slow = 1;
      kick(clevel, s2, s1, dt / 2, fast_or_slow, first_or_second, stime, etime);
    }

  fast_or_slow = 0;
  kick(clevel, s1, s2, dt / 2, fast_or_slow, first_or_second, stime, etime);

  drift(clevel, s1, etime, dt);
  first_or_second = 1;

  fast_or_slow = 0;
  kick(clevel, s1, s2, dt / 2, fast_or_slow, first_or_second, stime, etime);

  if(s2.n > 0)
    {
      fast_or_slow = 1;
      kick(clevel, s2, s1, dt / 2, fast_or_slow, first_or_second, stime, etime);
    }
}


void dkd(int clevel, struct sys s1, struct sys s2, real_t stime, real_t etime, real_t dt)
{
  int fast_or_slow = 0;
  int first_or_second = 1;

  drift(clevel, s1, etime - dt / 2, dt / 2);

  kick(clevel, s1, s2, dt, fast_or_slow, first_or_second, stime, etime);

  if(s2.n > 0)
    {
      fast_or_slow = 1;
      kick(clevel, s2, s1, dt, fast_or_slow, first_or_second, stime, etime);
    }

  drift(clevel, s1, etime, dt / 2);
}


static void split(real_t dt, struct sys s, struct sys *slow, struct sys *fast)
{
  real_t max_dt = -HUGE;
  real_t min_dt =  HUGE;

  size_t i = 0;
  struct particle *left, *right;

  left = s.part;
  right = s.last;
  //dt = fabs(dt);

  if(s.min_dt < dt) { //worth splitting if min_dt is less   than pivot time

    if(s.max_dt > dt) { //worth splitting if max_dt is larger than pivot time

#pragma omp parallel for reduction(max:max_dt) reduction(min:min_dt)
      for(size_t b=0; b<s.n; b++) {
	if(s.part[b].timestep > max_dt)
	  max_dt = s.part[b].timestep;

	if(s.part[b].timestep < min_dt)
	  min_dt = s.part[b].timestep;
      }

      while(1)
	{
	  if(i >= s.n)
	    printf("split error 1\n");
	  i++;
	  while(left->timestep < dt && left < right)
	    left++;
	  while(right->timestep >= dt && left < right)
	    right--;
	  if(left < right)
	    {
	      SWAP(*left, *right, struct particle);
	    }
	  else
	    break;
	}

      if(left->timestep < dt)
	left++;

      slow->n = s.last - left + 1;
      fast->n = (left - s.part);

      if(fast->n == 1)
	{
	  fast->n = 0;
	  slow->n = s.n;
	}
    }
    else {
      fast->n=s.n; //if s.max_dt < dt then fast is the entire subsystem
      slow->n=0;    
    }
  }
  else {
    fast->n=0; //if s.min_dt > dt then slow consistes of the entire subsystem
    slow->n=s.n;
  }

  if(slow->n > 0)
    {
      slow->part   = s.part + fast->n;
      slow->last   = s.last;
      slow->max_dt = s.max_dt;
      slow->min_dt = dt;
    }

  if(fast->n > 0)
    {
      fast->part = s.part;
      fast->last = s.part + fast->n - 1;
      if(s.max_dt > dt)
	fast->max_dt = max_dt;  //split is true, use new max_dt
      else
        fast->max_dt = s.max_dt;//no split, just use old max_dt
    }

  if(fast->n + slow->n != s.n)
    printf("split error 2, fast=%d slow=%d \n", fast->n, slow->n);
}

void
kick(int clevel, struct sys s1, struct sys s2, real_t dt, int fast_or_slow,
     int first_or_second, real_t t1, real_t t2)
{
  real_t boxsize = diag->BoxSize;
  real_t thistime = 1, a2inv = 1, a3inv = 1;

  thistime = exp((t1 + t2) / 2);

  a2inv = 1.0 / (thistime * thistime);
  a3inv = a2inv / thistime;

  real_t hubble_a = diag->Omega0 * a3inv + diag->OmegaLambda +
    (1.0 - diag->Omega0 - diag->OmegaLambda) * a2inv;

  hubble_a = Hubble * sqrt(hubble_a);

  real_t kick_factor = get_gravkick_factor(t1, t2) * G;

  //bodies.clear();
  bodies.resize(s1.n + s2.n);

  for(size_t i = 0; i < s1.n; i++)
    {
      for(int d = 0; d < 3; d++)
	bodies[i].X[d] = s1.part[i].pos[d];

      bodies[i].index = i;
      //bodies[i].q      = s1.part[i].mass;
      //bodies[i].type   = s1.part[i].type;
      if(fast_or_slow)
	bodies[i].issource = 0;
      else
	bodies[i].issource = 1;
      bodies[i].issink = 1;
      bodies[i].acc_old = s1.part[i].acc_mag;
    }

  for(size_t i = s1.n; i < s1.n + s2.n; i++)
    {
      for(int d = 0; d < 3; d++)
	bodies[i].X[d] = s2.part[i - s1.n].pos[d];

      bodies[i].index = i;
      //bodies[i].q        = s2.part[i-s1.n].mass;
      //bodies[i].type     = s2.part[i-s1.n].type;
      bodies[i].issource = 1;
      if(fast_or_slow || s1.n + s2.n < totnumBodies)
	bodies[i].issink = 0;
      else
	bodies[i].issink = 1;
      bodies[i].acc_old = s2.part[i - s1.n].acc_mag;
    }

  get_force_and_potential(bodies, 1, 1);

  float * force = new float[3 * bodies.size()];

  for(size_t b = 0; b < bodies.size(); b++)
    {
      unsigned int i = bodies[b].index;

      force[3 * i + 0] = bodies[b].d.F[0];
      force[3 * i + 1] = bodies[b].d.F[1];
      force[3 * i + 2] = bodies[b].d.F[2];
    }

  if(!fast_or_slow)
    {
      if(s1.n + s2.n == totnumBodies)
	{
	  real_t acc_mag;

#pragma omp parallel for
	  for(size_t i = 0; i < s1.n; i++)
	    {
	      acc_mag = sqrt(force[i * 3 + 0] * force[i * 3 + 0] +
			     force[i * 3 + 1] * force[i * 3 + 1] + 
			     force[i * 3 + 2] * force[i * 3 + 2]);

	      s1.part[i].acc_mag = acc_mag;

	      s1.part[i].timestep_all = sqrt(eta * thistime * softtype[1] / (acc_mag * a2inv * G));

	      s1.part[i].timestep_all *= hubble_a;

	      if(s1.part[i].timestep_all > MaxSizeTimestep)
		s1.part[i].timestep_all = MaxSizeTimestep;

	      s1.part[i].timestep = s1.part[i].timestep_all;
	    }

#pragma omp parallel for
	  for(size_t i = s1.n; i < s1.n + s2.n; i++)
	    {
	      acc_mag = sqrt(force[i * 3 + 0] * force[i * 3 + 0] +
			     force[i * 3 + 1] * force[i * 3 + 1] + 
			     force[i * 3 + 2] * force[i * 3 + 2]);

	      s2.part[i - s1.n].acc_mag = acc_mag;

	      s2.part[i - s1.n].timestep_all = sqrt(eta * thistime * softtype[1] / (acc_mag * a2inv * G));

	      s2.part[i - s1.n].timestep_all *= hubble_a;

	      if(s2.part[i - s1.n].timestep_all > MaxSizeTimestep)
		s2.part[i - s1.n].timestep_all = MaxSizeTimestep;

	      s2.part[i - s1.n].timestep = s2.part[i - s1.n].timestep_all;

	    }
	}
    }

  if(first_or_second)
    {
      if(s1.n + s2.n < totnumBodies)
	{
	  if(!fast_or_slow)
	    {
	      for(size_t i = 0; i < s1.n; i++)
		{
		  s1.part[i].acc_part[0] = force[i * 3 + 0];
		  s1.part[i].acc_part[1] = force[i * 3 + 1];
		  s1.part[i].acc_part[2] = force[i * 3 + 2];
		}
	    }
	  else
	    {
	      for(size_t i = 0; i < s1.n; i++)
		{
		  s1.part[i].acc_part[0] += force[i * 3 + 0];
		  s1.part[i].acc_part[1] += force[i * 3 + 1];
		  s1.part[i].acc_part[2] += force[i * 3 + 2];
		}
	    }

#pragma omp parallel for
	  for(size_t i = 0; i < s1.n; i++)
	    {
	      real_t acc_mag;

	      acc_mag =
		sqrt(s1.part[i].acc_part[0] * s1.part[i].acc_part[0] +
		     s1.part[i].acc_part[1] * s1.part[i].acc_part[1] +
		     s1.part[i].acc_part[2] * s1.part[i].acc_part[2]);

	      real_t newstep = sqrt(eta * thistime * softtype[1] / (acc_mag * a2inv * G));

	      newstep *= hubble_a;

	      if(newstep < s1.part[i].timestep_all)
		s1.part[i].timestep = newstep;
	      else
		s1.part[i].timestep = s1.part[i].timestep_all;
	    }
	}
    }

#pragma omp parallel for
  for(size_t i = 0; i < s1.n; i++)
    {
      s1.part[i].vel[0] += force[i * 3 + 0] * kick_factor;
      s1.part[i].vel[1] += force[i * 3 + 1] * kick_factor;
      s1.part[i].vel[2] += force[i * 3 + 2] * kick_factor;
    }

  delete[] force;
}

void drift(int clevel, struct sys s, real_t etime, real_t dt)
{
  real_t boxsize = diag->BoxSize;

  real_t drift_factor = (real_t) get_drift_factor(etime - dt, etime);

#pragma omp parallel for
  for(size_t i = 0; i < s.n; i++)
    {
      s.part[i].pos[0] += drift_factor * s.part[i].vel[0];
      s.part[i].pos[1] += drift_factor * s.part[i].vel[1];
      s.part[i].pos[2] += drift_factor * s.part[i].vel[2];

      if(s.part[i].pos[0] < 0)
	s.part[i].pos[0] += boxsize;
      if(s.part[i].pos[1] < 0)
	s.part[i].pos[1] += boxsize;
      if(s.part[i].pos[2] < 0)
	s.part[i].pos[2] += boxsize;

      if(s.part[i].pos[0] > boxsize)
	s.part[i].pos[0] -= boxsize;
      if(s.part[i].pos[1] > boxsize)
	s.part[i].pos[1] -= boxsize;
      if(s.part[i].pos[2] > boxsize)
	s.part[i].pos[2] -= boxsize;
    }
}


void evolve_split_hold(int clevel, struct sys s, real_t stime, real_t etime, real_t dt, int calc_timestep)
{
  struct sys slow = zerosys, fast = zerosys;

  split((real_t) dt, s, &slow, &fast);

  if(fast.n == 0)
    {
      diag->simtime += dt;
      real_t kick_factor = (real_t) get_gravkick_factor(diag->simtime - dt,
							diag->simtime);
      real_t drift_factor = (real_t) get_drift_factor(diag->simtime - dt,
						      diag->simtime);
      printf
	("---log(a)=%6.4f a=%6.4f redshift=%6.4f kick=%g drift=%g \n\n",
	 diag->simtime, exp(diag->simtime), 1.0/exp(diag->simtime)-1, 
	 kick_factor, drift_factor);
      fflush(stdout);
    }

  if(fast.n > 0)
    evolve_split_hold(clevel + 1, fast, stime, stime + dt / 2, dt / 2, 0);
  if(slow.n > 0)
    kdk(clevel, slow, fast, stime, etime, dt);
  if(fast.n > 0)
    evolve_split_hold(clevel + 1, fast, stime + dt / 2, etime, dt / 2, 1);
}


void evolve_split_hold_dkd(int clevel, struct sys s, real_t stime, real_t etime, real_t dt, int calc_timestep)
{
  struct sys slow = zerosys, fast = zerosys;

  start("split fast/slow");
  split((real_t) dt, s, &slow, &fast);
  stop("split fast/slow");
  
  if(fast.n == 0)
    {
      diag->simtime += dt;
      real_t kick_factor = (real_t) get_gravkick_factor(diag->simtime - dt,
							diag->simtime);
      real_t drift_factor = (real_t) get_drift_factor(diag->simtime - dt,
						      diag->simtime);

      printf
	("---log(a)=%6.4f a=%6.4f redshift=%6.4f \n kick=%g drift=%g \n\n",
	 diag->simtime, exp(diag->simtime), 1.0 / exp(diag->simtime)-1, 
	 kick_factor, drift_factor);
      fflush(stdout);
    }

  if(fast.n > 0)
    evolve_split_hold_dkd(clevel + 1, fast, stime, stime + dt / 2, dt / 2, 0);
  if(slow.n > 0)
    dkd(clevel, slow, fast, stime, etime, dt);
  if(fast.n > 0)
    evolve_split_hold_dkd(clevel + 1, fast, stime + dt / 2, etime, dt / 2, 1);
}

void evolve_individual(int clevel, struct sys s, real_t stime, real_t etime, real_t dt, int calc_timestep)
{
#if INDI
  real_t boxsize = diag->BoxSize;

  int TIMEBASE = 1 << 28;
  real_t basedt = dt / TIMEBASE;

#pragma omp parallel for
  for(size_t i = 0; i < s.n; i++)
    {
      s.part[i].Ti_begstep = 0;
      s.part[i].Ti_endstep = 0;
      s.part[i].timestep = 0;
    }

  int Ti_Current = 0;

  do
    {
      real_t t_drift;

      int min_glob = s.part[0].Ti_endstep;

      /* get the active number of particles for this timestep */
      for(size_t i = 1; i < s.n; i++)
	if(min_glob > s.part[i].Ti_endstep)
	  min_glob = s.part[i].Ti_endstep;

      diag->simtime += (min_glob - Ti_Current) * basedt;

      t_drift = (min_glob - Ti_Current) * basedt;

      real_t drift_factor = get_drift_factor(stime + Ti_Current * basedt,
					     stime + min_glob * basedt);

      Ti_Current = min_glob;

      real_t a2inv, a3inv, thistime;

      thistime = min_glob * basedt + stime;
      thistime = exp(thistime);

      a2inv = 1.0 / (thistime * thistime);
      a3inv = a2inv / thistime;

      real_t hubble_a = diag->Omega0 * a3inv + diag->OmegaLambda +
	(1.0 - diag->Omega0 - diag->OmegaLambda) * a2inv;

      hubble_a = Hubble * sqrt(hubble_a);

      printf("---log(a)=%6.4f a=%6.4f redshift=%6.4f \n\n", 
	     log(thistime), thistime, 1.0 / thistime-1);
      fflush(stdout);

      int activecount = 0;

#pragma omp parallel for
      for(size_t i = 0; i < s.n; i++)
	{
	  if(s.part[i].Ti_endstep == Ti_Current)
	    {
	      s.part[i].isactive = 1;
	      activecount++;
	    }
	  else
	    s.part[i].isactive = 0;

	  s.part[i].pos[0] += drift_factor * s.part[i].vel[0];
	  s.part[i].pos[1] += drift_factor * s.part[i].vel[1];
	  s.part[i].pos[2] += drift_factor * s.part[i].vel[2];

	  if(s.part[i].pos[0] < 0)
	    s.part[i].pos[0] += boxsize;
	  if(s.part[i].pos[1] < 0)
	    s.part[i].pos[1] += boxsize;
	  if(s.part[i].pos[2] < 0)
	    s.part[i].pos[2] += boxsize;

	  if(s.part[i].pos[0] > boxsize)
	    s.part[i].pos[0] -= boxsize;
	  if(s.part[i].pos[1] > boxsize)
	    s.part[i].pos[1] -= boxsize;
	  if(s.part[i].pos[2] > boxsize)
	    s.part[i].pos[2] -= boxsize;
	}

      //bodies.clear();
      bodies.resize(s.n);

      for(size_t i = 0; i < s.n; i++)
	{
	  for(int d = 0; d < 3; d++)
	    {
	      bodies[i].X[d] = s.part[i].pos[d];
	    }
	  if(s.part[i].isactive < 1)
	    bodies[i].issink = 0;
	  else
	    bodies[i].issink = 1;
	  bodies[i].index = i;
	  //bodies[i].type    = s.part[i].type;
	  //bodies[i].q       = s.part[i].mass;
	  bodies[i].issource = 1;
	  bodies[i].acc_old = s.part[i].acc_mag;
	}

      get_force_and_potential(bodies, 0, 1);

      int ti_min;

#pragma omp parallel for
      for(size_t b = 0; b < bodies.size(); b++)
	{
	  unsigned int i = bodies[b].index;

	  if(s.part[i].isactive < 1)
	    continue;

	  real_t acc_mag = sqrt(bodies[b].d.F[0] * bodies[b].d.F[0] +
				bodies[b].d.F[1] * bodies[b].d.F[1] + bodies[b].d.F[2] * bodies[b].d.F[2]);

	  s.part[i].acc_mag = acc_mag;

	  real_t timestep = sqrt(eta * thistime * softtype[1] / (acc_mag * a2inv * G));

	  timestep *= hubble_a;

	  if(timestep > MaxSizeTimestep)
	    timestep = MaxSizeTimestep;

	  int ti_step = timestep / basedt;

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

	  if(Ti_Current > 0)
	    {
	      if(ti_step > (s.part[i].Ti_endstep - s.part[i].Ti_begstep))	/* timestep wants to increase */
		{
		  if(((TIMEBASE - s.part[i].Ti_endstep) % ti_step) > 0)
		    ti_step = s.part[i].Ti_endstep - s.part[i].Ti_begstep;	/* leave at old step */
		}
	    }

	  if(Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	    ti_step = 0;

	  if((TIMEBASE - Ti_Current) < ti_step)
	    ti_step = TIMEBASE - Ti_Current;

	  int tstart = (s.part[i].Ti_begstep + s.part[i].Ti_endstep) / 2;
	  int tend = s.part[i].Ti_endstep + ti_step / 2;

	  s.part[i].timestep = ti_step * basedt;

	  s.part[i].Ti_begstep = s.part[i].Ti_endstep;
	  s.part[i].Ti_endstep = s.part[i].Ti_begstep + ti_step;

	  real_t kick_factor = get_gravkick_factor(stime + tstart * basedt,
						   stime + tend * basedt) * G;

	  s.part[i].vel[0] += bodies[b].d.F[0] * kick_factor;
	  s.part[i].vel[1] += bodies[b].d.F[1] * kick_factor;
	  s.part[i].vel[2] += bodies[b].d.F[2] * kick_factor;
	}
    }
  while(Ti_Current < TIMEBASE);
#endif
}

void do_evolve(struct sys s, double dt)
{
  int clevel;

  if(dt == 0)
    return;
  clevel = 0;

#if INDI
  evolve_individual(clevel, s, (real_t) diag->simtime, (real_t) diag->simtime + dt, (real_t) dt, 1);
#endif

#if HOLD
  evolve_split_hold_dkd(clevel, s, (real_t) diag->simtime, (real_t) diag->simtime + dt, (real_t) dt, 1);
#endif

  //evolve_split_hold(clevel,s,(real_t)diag->simtime,(real_t)diag->simtime+dt, (real_t) dt, 1);
}

#include "cosmo.h"

#endif
