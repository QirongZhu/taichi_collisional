#ifndef integrator_h
#define integrator_h

#include "kernel_fast_lean.h"
#include "traverse_eager.h"
#include "timer.h"

#define MAXLEVEL 64

using namespace exafmm;

extern struct sys zerosys;

/* diagnostics */
struct diagnostics
{
  double simtime;
};

extern struct diagnostics global_diag;
extern struct diagnostics *diag;
#pragma omp threadprivate(diag)

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}
#define ABS(X) (((X) >= 0) ? (X) : -(X))
#define SIGN(X)   ((X>0)-(X<0))
#define LOG(fmt, ...) {				\
    printf("%s:%d\t", __FILE__, __LINE__);	\
    printf(fmt, ## __VA_ARGS__);		\
  }


#define CHECK_TIMESTEP(etime,stime,dt,clevel)				\
  if(sizeof(dt)==sizeof(long double)) {					\
    if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)			\
      ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel); \
  } else {								\
    if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)			\
      ENDRUN("timestep too small: etime=%le stime=%le dt=%le clevel=%u\n", (double) etime, (double) stime, (double) dt, clevel); \
  }

#define COMPSUM(s,e,ds){double a=s;e+=ds;s=(a+e);e+=(a-s);}

void get_force_and_potential(Bodies & bodies, bool get_steps);

void split(double dt, struct sys s, struct sys *slow, struct sys *fast);

void kick_cpu(int clevel, struct sys s1, struct sys s2, double dt);

void kick(int clevel, struct sys sinks, struct sys sources, double dt,
	  bool update_timestep, bool sinks_are_fast);
/* =kick sys1 for interactions with sys2  */

void drift(int clevel, struct sys s, double etime, double dt);	/* drift sys */

void drift_naive(int clevel, struct sys s, double etime);

void evolve_split_naive(int clevel, struct sys sys1, struct sys sys2,
			double stime, double etime, double dt, int calc_timestep);

void evolve_split_hold_dkd(int clevel, struct sys s, double stime,
			   double etime, double dt, int calc_timestep);


void kick_cpu(int clevel, struct sys s1, struct sys s2, double dt)
{
  unsigned int i, j;
  real_t dx[3], dr3, dr2, dr, acci;
  real_t acc[3];

  for(i = 0; i < s1.n; i++)
    {
      acc[0] = 0.;
      acc[1] = 0.;
      acc[2] = 0.;

      for(j = 0; j < s2.n; j++)
	{
	  dx[0] = s1.part[i].pos[0] - s2.part[j].pos[0];
	  dx[1] = s1.part[i].pos[1] - s2.part[j].pos[1];
	  dx[2] = s1.part[i].pos[2] - s2.part[j].pos[2];
	  dr2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
	  if(dr2 > 0)
	    {
	      dr = sqrt(dr2);
	      dr3 = dr * dr2;
	      acci = s2.part[j].mass / dr3;

	      acc[0] -= dx[0] * acci;
	      acc[1] -= dx[1] * acci;
	      acc[2] -= dx[2] * acci;
	    }
	}
      s1.part[i].vel[0] += dt * acc[0];
      s1.part[i].vel[1] += dt * acc[1];
      s1.part[i].vel[2] += dt * acc[2];
    }
}

void kick_naive(int rung, struct sys sinks, struct sys sources1, struct sys sources2,
		double dt, bool update_timestep)
{

#if DEBUG
  start("Prepare");
#endif

  bodies.clear();
  bodies.resize(sinks.n + sources1.n + sources2.n);

  for(unsigned int i = 0; i < sinks.n; i++)
    {
      for(int d = 0; d < 3; d++)
	{
	  bodies[i].X[d] = sinks.part[i].pos[d];
	  bodies[i].V[d] = sinks.part[i].vel[d];
	}
      bodies[i].index = i;
      bodies[i].q = sinks.part[i].mass;
      bodies[i].issink = true;
      bodies[i].issource = true;
      bodies[i].timestep = 0;
    }

  if(sources1.n > 0)
    {
      for(unsigned int i = sinks.n; i < sinks.n + sources1.n; i++)
	{
	  for(int d = 0; d < 3; d++)
	    {
	      bodies[i].X[d] = sources1.part[i - sinks.n].pos[d];
	      bodies[i].V[d] = sources1.part[i - sinks.n].vel[d];
	    }
	  bodies[i].index = i;
	  bodies[i].q = sources1.part[i - sinks.n].mass;
	  bodies[i].issource = true;
	  bodies[i].issink = false;
	}
    }

  if(sources2.n > 0)
    {
      unsigned int startindex = sinks.n + sources1.n;
      for(unsigned int i = startindex; i < startindex + sources2.n; i++)
	{
	  for(int d = 0; d < 3; d++)
	    {
	      bodies[i].X[d] = sources2.part[i - startindex].pos[d];
	      bodies[i].V[d] = sources2.part[i - startindex].vel[d];
	    }

	  bodies[i].index = i;
	  bodies[i].q = sources2.part[i - startindex].mass;
	  bodies[i].issource = true;
	  bodies[i].issink = false;
	}
    }

#if DEBUG
  stop("Prepare");
#endif

  get_force_and_potential(bodies, update_timestep);

#if DEBUG
  start("Getting result");
#endif

  real_t *force = new real_t[3 * bodies.size()];
  real_t *potential = new real_t[bodies.size()];
  real_t *timestep = new real_t[bodies.size()];
  real_t *acc_old = new real_t[bodies.size()];

  for(size_t b = 0; b < bodies.size(); b++)
    {
      unsigned int i = bodies[b].index;
      force[3 * i + 0] = bodies[b].F[0];
      force[3 * i + 1] = bodies[b].F[1];
      force[3 * i + 2] = bodies[b].F[2];
      potential[i] = bodies[b].p;
      timestep[i]  = (real_t)1/sqrt(sqrt(bodies[b].timestep));
      acc_old[i]   = bodies[b].acc_old;
    }

  if(update_timestep)
    {
      for(size_t i = 0; i < sinks.n; i++)
	{
	  sinks.part[i].timestep = timestep[i];
	}
    }

  for(size_t i = 0; i < sinks.n; i++)
    {
      sinks.part[i].pot = potential[i];
      sinks.part[i].vel[0] += force[i * 3 + 0] * dt;
      sinks.part[i].vel[1] += force[i * 3 + 1] * dt;
      sinks.part[i].vel[2] += force[i * 3 + 2] * dt;

      sinks.part[i].acc_old = acc_old[i];

#ifdef OUTPUTACC
      sinks.part[i].acc[0] = force[i * 3 + 0];
      sinks.part[i].acc[1] = force[i * 3 + 1];
      sinks.part[i].acc[2] = force[i * 3 + 2];
#endif
    }

  delete[]force;
  delete[]potential;
  delete[]timestep;

#if DEBUG
  stop("Getting result");
#endif
}

void kick(int clevel, struct sys sinks, struct sys sources,
	  double dt, bool update_timestep, bool sinks_are_fast)
{

  bodies.resize(sinks.n + sources.n);

  for(size_t i = 0; i < sinks.n; i++)
    {
      for(int d = 0; d < 3; d++)
	{
	  bodies[i].X[d] = sinks.part[i].pos[d];
	  //bodies[i].V[d] = sinks.part[i].vel[d];
	}

      bodies[i].index = i;
      bodies[i].q = sinks.part[i].mass;
      if(sinks_are_fast)
	bodies[i].issource = false;
      else
	bodies[i].issource = true;
      bodies[i].issink = true;
      bodies[i].p = 0;
      bodies[i].F[0] = 0;
      bodies[i].F[1] = 0;
      bodies[i].F[2] = 0;
      //bodies[i].timestep   = 0;
    }

  for(size_t i = sinks.n; i < sinks.n + sources.n; i++)
    {
      for(int d = 0; d < 3; d++)
	{
	  bodies[i].X[d] = sources.part[i - sinks.n].pos[d];
	  //bodies[i].V[d]   = sources.part[i-sinks.n].vel[d];
	}
      bodies[i].index = i;
      bodies[i].q = sources.part[i - sinks.n].mass;
      bodies[i].issource = true;
      bodies[i].issink = false;
      bodies[i].p = 0;
      bodies[i].F[0] = 0;
      bodies[i].F[1] = 0;
      bodies[i].F[2] = 0;
      //  bodies[i].timestep = 0;
    }

  get_force_and_potential(bodies, false);

  real_t *force = new real_t[3 * bodies.size()];
  real_t *potential = new real_t[bodies.size()];
  //real_t * timestep  = new real_t[  bodies.size()];

  for(size_t b = 0; b < bodies.size(); b++)
    {
      unsigned int i = bodies[b].index;
      force[3 * i + 0] = bodies[b].F[0];
      force[3 * i + 1] = bodies[b].F[1];
      force[3 * i + 2] = bodies[b].F[2];
      potential[i] = bodies[b].p;
      //timestep[i]  = bodies[b].timestep;
    }

  // if(update_timestep)
  //   for(size_t i = 0; i < sinks.n; i++)
  //   { 
  //     sinks.part[i].timestep = timestep[i];
  //   }

  for(size_t i = 0; i < sinks.n; i++)
    {
      sinks.part[i].pot = potential[i];
      sinks.part[i].vel[0] += force[i * 3 + 0] * dt;
      sinks.part[i].vel[1] += force[i * 3 + 1] * dt;
      sinks.part[i].vel[2] += force[i * 3 + 2] * dt;
    }

  delete[]force;
  delete[]potential;
  // delete[] timestep;
}

void get_force_and_potential(Bodies & bodies, bool get_steps)
{
#ifdef FMM
  if(bodies.size() < 8e3)
    {
      direct(bodies, bodies, get_steps);
    }
  else
    {
#if DEBUG
      start("buildTree");
#endif
      Cells cells = buildTree(bodies);

#if DEBUG
      stop("buildTree");
#endif

#if DEBUG
      start("upwardPass_low");
#endif
      upwardPass_low(cells);
#if DEBUG
      stop("upwardPass_low");
#endif
#if DEBUG
      start("horizontalPass_low");
#endif
      horizontalPass_low(cells, cells);
#if DEBUG
      stop("horizontalPass_low");
#endif

#if DEBUG
      start("downwardPass");
#endif
      downwardPass_low(cells);
#if DEBUG
      stop("downwardPass");
#endif

      for(size_t b = 0; b < bodies.size(); b++)
	{
	  bodies[b].F[0] = 0;
	  bodies[b].F[1] = 0;
	  bodies[b].F[2] = 0;
	  bodies[b].p = 0;
	  bodies[b].timestep = 0;
	}

      for(size_t i = 0; i < cells.size(); i++)
	{
	  for(int term = 0; term < NTERM; term++)
	    {
	      cells[i].L[term] = 0;
	      cells[i].M[term] = 0;
	      cells[i].min_acc = HUGE;
	      cells[i].has_sink = false;
	    }
	}

#if DEBUG
      start("upwardPass");
#endif
      upwardPass(cells);
#if DEBUG
      stop("upwardPass");
#endif
      for(size_t i = 0; i < cells.size(); i++)
	{
	  for(int n = 0; n <= P; n++)
	    {
	      for(int m = -n; m <= n; m++)
		{
		  int index = n * n + n + m;
		  cells[i].Pn[n] += cells[i].M[index] *
		    cells[i].M[index] * factorial_table[n + m] * factorial_table[n - m];
		}
	      cells[i].Pn[n] = sqrt(cells[i].Pn[n]);
	    }
	}
#if DEBUG
      start("horizontalPass");
#endif
      horizontalPass(cells, cells, 1, get_steps);
#if DEBUG
      stop("horizontalPass");
#endif
#if DEBUG
      start("downwardPass");
#endif
      downwardPass(cells);
#if DEBUG
      stop("downwardPass");
#endif

    }

#else // if FMM not defined, use direct summation instead

  direct(bodies, bodies, get_steps);

#endif
}

struct sys join(struct sys sinks, struct sys sources)
{
  struct sys s = zerosys;
  if(sinks.n == 0)
    return sources;
  if(sources.n == 0)
    return sinks;
  s.n = sinks.n + sources.n;
  if(sinks.part + sinks.n == sources.part)
    {
      s.part = sinks.part;
      s.last = sources.last;
    }
  else
    {
      if(sources.part + sources.n == sinks.part)
	{
	  s.part = sources.part;
	  s.last = sinks.last;
	}
      else
	ENDRUN("join error 1");
    }
  if(s.last - s.part + 1 != s.n)
    ENDRUN("join error 2");
  return s;
}

void split(double dt, struct sys s, struct sys *slow, struct sys *fast)
{
  unsigned int i = 0;
  struct particle *left, *right;
  left = s.part;
  right = s.last;
  dt = fabs(dt);
  while(1)
    {
      if(i >= s.n)
	ENDRUN("split error 1");
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

  if(slow->n > 0)
    {
      slow->part = s.part + fast->n;
      slow->last = s.last;
    }

  if(fast->n > 0)
    {
      fast->part = s.part;
      fast->last = s.part + fast->n - 1;
    }

  if(fast->n + slow->n != s.n)
    ENDRUN("split error 2");
}

void drift(int clevel, struct sys s, double etime, double dt)
{
  for(unsigned int i = 0; i < s.n; i++)
    {
      COMPSUM(s.part[i].pos[0], s.part[i].pos_e[0], dt * s.part[i].vel[0]);
      COMPSUM(s.part[i].pos[1], s.part[i].pos_e[1], dt * s.part[i].vel[1]);
      COMPSUM(s.part[i].pos[2], s.part[i].pos_e[2], dt * s.part[i].vel[2]);
      //  s.part[i].pos[0] += dt*s.part[i].vel[0];
      //  s.part[i].pos[1] += dt*s.part[i].vel[1];
      //  s.part[i].pos[2] += dt*s.part[i].vel[2];
      s.part[i].postime = etime;
    }
}

void drift_naive(int clevel, struct sys s, double etime)
{
  double dt;
  for(unsigned int i = 0; i < s.n; i++)
    {
      dt = etime - s.part[i].postime;
      COMPSUM(s.part[i].pos[0], s.part[i].pos_e[0], dt * s.part[i].vel[0]);
      COMPSUM(s.part[i].pos[1], s.part[i].pos_e[1], dt * s.part[i].vel[1]);
      COMPSUM(s.part[i].pos[2], s.part[i].pos_e[2], dt * s.part[i].vel[2]);
      //s.part[i].pos[0] += dt*s.part[i].vel[0];
      //s.part[i].pos[1] += dt*s.part[i].vel[1];
      //s.part[i].pos[2] += dt*s.part[i].vel[2];
      s.part[i].postime = etime;
    }
}

void evolve_split_naive(int clevel, struct sys sys1, struct sys sys2,
			double stime, double etime, double dt, int calc_timestep)
{
  struct sys slow = zerosys, fast = zerosys;

  //CHECK_TIMESTEP(etime,stime,dt,clevel);

  split(dt, sys1, &slow, &fast);

  if(fast.n == 0)
    {
      diag->simtime += dt;
    }

  bool update_timestep = false;
  bool sinks_are_fast = false;

  if(slow.n > 0)
    kick_naive(clevel, slow, fast, sys2, dt / 2, update_timestep);

  //kick(clevel, slow, join(sys1,sys2), dt/2, update_timestep, sinks_are_fast);

  if(fast.n > 0)
    evolve_split_naive(clevel + 1, fast, join(slow, sys2), stime, stime + dt / 2, dt / 2, 0);

  if(slow.n > 0)
    drift_naive(clevel, join(slow, sys2), stime + dt / 2);

  if(fast.n > 0)
    evolve_split_naive(clevel + 1, fast, join(slow, sys2), stime + dt / 2, etime, dt / 2, 1);

  if(slow.n > 0)
    drift_naive(clevel, join(slow, sys2), etime);

  update_timestep = true;
  if(slow.n > 0)
    if(slow.n > 0)
      kick_naive(clevel, slow, fast, sys2, dt / 2, update_timestep);

  //kick(clevel, slow, join(sys1,sys2), dt/2, update_timestep, sinks_are_fast);
}

void evolve_split_hold_dkd(int clevel, struct sys s, double stime, double etime, double dt,
			   bool calc_timestep)
{
  struct sys slow = zerosys, fast = zerosys;

  if(calc_timestep){    
    kick_naive(clevel, s, zerosys, zerosys, 0, true);
  }

  split(dt, s, &slow, &fast);

  if(fast.n == 0)
    {
      diag->simtime += dt;
#if DEBUG
      printf("t=%g s=%d \n", diag->simtime, s.n);
#endif
    }

  fflush(stdout);

  if(fast.n > 0)
    evolve_split_hold_dkd(clevel + 1, fast, stime, stime + dt / 2, dt / 2, false);

  drift(clevel, slow, stime + dt / 2, dt / 2);

  bool update_timestep = false;

  if(s.n > 128)
    {
      bool sinks_are_fast = false;
      kick(clevel, slow, fast, dt, update_timestep, sinks_are_fast);
      
      sinks_are_fast = true;
      if(fast.n > 0)
	kick(clevel, fast, slow, dt, update_timestep, sinks_are_fast);
    }
  else if(slow.n > 0)
    {
      kick_cpu(clevel, slow, join(slow, fast), dt);
      if(fast.n > 0)
    	kick_cpu(clevel, fast, slow, dt);
    }

  drift(clevel, slow, etime, dt / 2);

  if(fast.n > 0)
    evolve_split_hold_dkd(clevel + 1, fast, stime + dt / 2, etime, dt / 2, true);
}

void do_evolve(struct sys s, double dt)
{
  int clevel;

  if(dt == 0)
    return;
  clevel = 0;

  kick_naive(clevel, s, zerosys, zerosys, 0, true);

#if BLOCK
  evolve_split_naive(clevel, s, zerosys, (real_t) diag->simtime, (real_t) diag->simtime + dt, (real_t) dt, 1);
#endif

#if HOLD
  evolve_split_hold_dkd(clevel, s, (real_t) diag->simtime, (real_t) diag->simtime + dt, (real_t) dt, true);
#endif
}

double system_kinetic_energy(struct sys s)
{
  unsigned int i;
  long double e = 0.;
  for(i = 0; i < s.n; i++)
    e += s.part[i].mass / G * (s.part[i].vel[0] * s.part[i].vel[0] +
			       s.part[i].vel[1] * s.part[i].vel[1] + 
			       s.part[i].vel[2] * s.part[i].vel[2]) / 2;
  return (double) e;
}

double system_potential_energy(struct sys s)
{
  unsigned int i;
  long double e = 0.;
  for(i = 0; i < s.n; i++)
    e += s.part[i].mass / G * s.part[i].pot / 2;
  return (double) e;
}

void system_center_of_mass(struct sys s, long double *cmpos, long double *cmvel)
{
  long double mass = 0., pos[3] = { 0., 0., 0. }, vel[3] =
  {
  0., 0., 0.};
  for(unsigned int p = 0; p < s.n; p++)
    {
      for(int i = 0; i < 3; i++)
	{
	  pos[i] += (long double) s.part[p].mass * s.part[p].pos[i];
	  vel[i] += (long double) s.part[p].mass * s.part[p].vel[i];
	}
      mass += (long double) s.part[p].mass;
    }
  for(int i = 0; i < 3; i++)
    {
      cmpos[i] = pos[i] / mass;
      cmvel[i] = vel[i] / mass;
    }
}

#endif
