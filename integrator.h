#pragma once

#include "kernel_fast_lean.h"
#include "traverse_eager.h"
#include "timestep.h"
#include "timer.h"
#include "universal.h"

#define MAXLEVEL 64

using namespace exafmm;

extern struct sys zerosys;

/* diagnostics */
struct diagnostics {
  double simtime;
};

extern struct diagnostics global_diag;
extern struct diagnostics *diag;
#pragma omp threadprivate(diag)

#define SWAP(a,b,c) {c t;t=(a);(a)=(b);(b)=t;}
#define ABS(X) (((X) >= 0) ? (X) :-(X))
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

#define COMPSUM(s,e,ds){real_t a=s;e+=ds;s=(a+e);e+=(a-s);}

void get_force_and_potential(Bodies & bodies);

void split(double dt, struct sys s, struct sys *slow, struct sys *fast);

void kick_cpu(int clevel, struct sys s1, struct sys s2, double dt,
              double b, double c, bool isgradient);

void kick(int clevel, struct sys sinks, struct sys sources, double dt,
	  bool update_timestep, bool sinks_are_fast);
/* =kick sys1 for interactions with sys2  */

void drift(int clevel, struct sys s, double etime, double dt);	/* drift sys */

void drift_naive(int clevel, struct sys s, double etime);

void evolve_split_naive(int clevel, struct sys sys1, struct sys sys2,
			double stime, double etime, double dt, int calc_timestep);

void evolve_split_hold_dkd(int clevel, struct sys s, double stime,
			   double etime, double dt, int calc_timestep);

void integrate_sym_six(int clevel, struct sys total,
                       double stime, double etime, double dt);

void integrate_sym_eight(int clevel, struct sys total,
			 double stime, double etime, double dt);

void evolve_kepler(int clevel, struct sys s, double stime, double etime, double dt) {

  if (s.n != 2) ENDRUN("two-body solver was called with sys.n=%u\n", s.n);

  // translate coordinates original frame to 2-body frame
  int k;

  double dpos[3],dpos0[3];
  double pos_cm[3];
  double dvel[3],dvel0[3];
  double vel_cm[3];
  double m1 = s.part->mass;
  double m2 = s.last->mass;
  double mtot = s.part->mass+s.last->mass;
  double f1 = m2 / mtot;
  double f2 = m1 / mtot;
  if(mtot>0.) {
      
    for(int d=0; d<3; d++){
      dpos0[d] = s.part->pos[d]-s.last->pos[d];
      dvel0[d] = s.part->vel[d]-s.last->vel[d];
    }
      
    for(int d=0; d<3; d++){
      pos_cm[d] = (m1*s.part->pos[d] + m2*s.last->pos[d]) / mtot;
      vel_cm[d] = (m1*s.part->vel[d] + m2*s.last->vel[d]) / mtot;
    }
    // evolve center of mass for dt
    for(int d=0; d<3; d++)
      pos_cm[d] += vel_cm[d] * dt;

    // call the Kepler solver
    Binary state0 = {dpos0[0],dpos0[1],dpos0[2], dvel0[0],dvel0[1],dvel0[2], 0, 0, 0};
    Binary state1;
              
    kepler_step(mtot, dt, &state0, &state1);
                  
    // translate coordinates from 2-body frame to original frame
    s.part->pos[0] = pos_cm[0]+f1*state1.x;
    s.part->pos[1] = pos_cm[1]+f1*state1.y;
    s.part->pos[2] = pos_cm[2]+f1*state1.z;
    s.part->vel[0] = vel_cm[0]+f1*state1.xd;
    s.part->vel[1] = vel_cm[1]+f1*state1.yd;
    s.part->vel[2] = vel_cm[2]+f1*state1.zd;
    s.last->pos[0] = pos_cm[0]-f2*state1.x;
    s.last->pos[1] = pos_cm[1]-f2*state1.y;
    s.last->pos[2] = pos_cm[2]-f2*state1.z;
    s.last->vel[0] = vel_cm[0]-f2*state1.xd;
    s.last->vel[1] = vel_cm[1]-f2*state1.yd;
    s.last->vel[2] = vel_cm[2]-f2*state1.zd;
      
  } else {
    for(int d=0; d<3; d++) {
      s.part->pos[d] +=s.part->vel[d] * dt;
      s.last->pos[d] +=s.last->vel[d] * dt;
    }
  }
}

void kick(int clevel, struct sys s, double etime, double dt)
{
#pragma omp parallel for if(s.n > ncrit)
  for(unsigned int i = 0; i < s.n; i++)
    {
      for(int d=0; d<3; d++) {
	COMPSUM(s.part[i].vel[d], s.part[i].vel_e[d], dt*s.part[i].acc[d]);
      }
      s.part[i].postime = etime;
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
      for(int d=0; d<3; d++){
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
	  for(int d=0; d<3; d++){
	    bodies[i].X[d] = sources1.part[i-sinks.n].pos[d];
	    bodies[i].V[d] = sources1.part[i-sinks.n].vel[d];
	  }
	  bodies[i].index = i;
	  bodies[i].q = sources1.part[i-sinks.n].mass;
	  bodies[i].issource = true;
	  bodies[i].issink = false;
	}
    }

  if(sources2.n > 0)
    {
      unsigned int startindex = sinks.n + sources1.n;
      for(unsigned int i = startindex; i < startindex + sources2.n; i++)
	{
	  for(int d=0; d<3; d++){
	    bodies[i].X[d] = sources2.part[i-startindex].pos[d];
	    bodies[i].V[d] = sources2.part[i-startindex].vel[d];
	  }
	  bodies[i].index = i;
	  bodies[i].q = sources2.part[i-startindex].mass;
	  bodies[i].issource = true;
	  bodies[i].issink = false;
	}
    }

#if DEBUG
  stop("Prepare");
#endif

  get_force_and_potential(bodies);

#if DEBUG
  start("Getting result");
#endif

  std::vector<real_t> force(bodies.size()*3);
  std::vector<real_t> potential(bodies.size());
  std::vector<real_t> acc_old(bodies.size());

  real_t min_timesteps = HUGE;
    
  for(size_t b = 0; b < bodies.size(); b++) {
    unsigned int i = bodies[b].index;
    force[3*i+0]   = bodies[b].F[0];
    force[3*i+1]   = bodies[b].F[1];
    force[3*i+2]   = bodies[b].F[2];
    potential[i] = bodies[b].p;
    acc_old[i]   = bodies[b].acc_old;
  }

  for(size_t i = 0; i < sinks.n; i++) {
    sinks.part[i].pot = potential[i];
    sinks.part[i].vel[0] += force[3*i+0]* dt;
    sinks.part[i].vel[1] += force[3*i+1]* dt;
    sinks.part[i].vel[2] += force[3*i+2]* dt;
    sinks.part[i].acc_old = acc_old[i];
  }

#if DEBUG
  stop("Getting result");
#endif
}

void kick_cpu(int clevel, struct sys s1, struct sys s2,
	      double dt, double b, double c, bool isgradient)
{

  real_t fac = 0;
    
  if(isgradient)
    fac = 2.0 * c / b * dt * dt;
    
#pragma omp parallel for if(s1.n > ncrit)
  for(unsigned int i = 0; i < s1.n; i++)
    {
      real_t acc[3] = {0.0, 0.0, 0.0};

      for(unsigned int j = 0; j < s2.n; j++)
	{
	  if(s1.part[i].id == s2.part[j].id)
	    continue;
        
	  real_t dr3, dr2, dr;

	  real_t dx[3];
	  for(int d=0; d<3; d++)
	    dx[d] = (s1.part[i].pos[d]-s2.part[j].pos[d])+(s1.part[i].acc[d]-s2.part[j].acc[d])*fac;
        
	  dr2  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	  dr   = sqrt(dr2);
	  dr3  = dr*dr2;
	  dr   = s2.part[j].mass / dr3;
        
	  for(int d=0; d<3; d++)
            acc[d] -= dx[d]*dr;
        
	}
        
      for(int d=0; d<3; d++){
	s1.part[i].acc[d]  = acc[d];
      }
        
      for(int d=0; d<3; d++) {
        COMPSUM(s1.part[i].vel[d], s1.part[i].vel_e[d], (b*dt)*s1.part[i].acc[d]);
      }

    }
}

void kick_self(int clevel, struct sys sinks, double dt,
	       double b, double c, bool isgradient)
{
  real_t fac = 0;

  if(isgradient)
    fac = 2.0 * c / b * dt * dt;

  bodies.resize(sinks.n);

  for(size_t i = 0; i < sinks.n; i++)
    {
      for(int d=0; d<3; d++){
	bodies[i].X[d]     = sinks.part[i].pos[d] + fac * sinks.part[i].acc[d];
      }
      bodies[i].index = i;
      bodies[i].q = sinks.part[i].mass;
      bodies[i].issource = true;
      bodies[i].issink = true;
      bodies[i].p    = 0;
      bodies[i].F[0] = 0;
      bodies[i].F[1] = 0;
      bodies[i].F[2] = 0;
    }

  get_force_and_potential(bodies);
        
  std::vector<real_t> force(bodies.size()*3);
  std::vector<real_t> potential(bodies.size());

  for(size_t b = 0; b < bodies.size(); b++)
    {
      unsigned int i = bodies[b].index;
      force[3*i+0]   = bodies[b].F[0];
      force[3*i+1]   = bodies[b].F[1];
      force[3*i+2]   = bodies[b].F[2];

      potential[i] = bodies[b].p;
    }
        
  if(dt > 0) {
    for(size_t i = 0; i < sinks.n; i++)
      {
	sinks.part[i].pot = potential[i];
	for(int d=0; d<3; d++) {
	  COMPSUM(sinks.part[i].vel[d], sinks.part[i].vel_e[d], (b*dt)*force[3*i+d]);
	}
      }
  }else{
    for(size_t i = 0; i < sinks.n; i++)
      {
	sinks.part[i].pot = potential[i];
	for(int d=0; d<3; d++) {
	  sinks.part[i].acc[d] = force[3 * i+d];
	}
      }
  }
    
}

void kick_sf(int clevel, struct sys sinks, struct sys sources,
	     double dt, double b, double c, bool isgradient)
{

  real_t fac =  0;
    
  if(isgradient)
    fac = 2.0 * c / b * dt * dt;

  bodies.resize(sinks.n + sources.n);

  for(size_t i = 0; i < sinks.n; i++)
    {
      for(int d=0; d<3; d++){
	bodies[i].X[d] = sinks.part[i].pos[d] + fac * sinks.part[i].acc[d];
      }
      bodies[i].index = i;
      bodies[i].q = sinks.part[i].mass;
      bodies[i].issource = false;
      bodies[i].issink = true;
      bodies[i].p = 0;
      bodies[i].F[0] = 0;
      bodies[i].F[1] = 0;
      bodies[i].F[2] = 0;
    }

  for(size_t i = sinks.n; i < sinks.n + sources.n; i++)
    {
      for(int d=0; d<3; d++){
	bodies[i].X[d] = sources.part[i-sinks.n].pos[d] + fac * sources.part[i-sinks.n].acc[d];
      }
      bodies[i].index = i;
      bodies[i].q = sources.part[i-sinks.n].mass;
      bodies[i].issource = true;
      bodies[i].issink = false;
      bodies[i].p = 0;
      bodies[i].F[0] = 0;
      bodies[i].F[1] = 0;
      bodies[i].F[2] = 0;
    }

  get_force_and_potential(bodies);
        
  std::vector<real_t> force(bodies.size()*3);
    
  std::vector<real_t> potential(bodies.size());

  for(size_t b = 0; b < bodies.size(); b++)
    {
      unsigned int i = bodies[b].index;
      force[3*i+0]   = bodies[b].F[0];
      force[3*i+1]   = bodies[b].F[1];
      force[3*i+2]   = bodies[b].F[2];

      potential[i] = bodies[b].p;
    }
        
  if(b > 0) {
    for(size_t i = 0; i < sinks.n; i++)
      {
	sinks.part[i].pot = potential[i];
	for(int d=0; d<3; d++)
	  COMPSUM(sinks.part[i].vel[d], sinks.part[i].vel_e[d], (b*dt) * force[3*i+d] );
      }
  }else{
    for(size_t i = 0; i < sinks.n; i++)
      {
	sinks.part[i].pot = potential[i];
	for(int d=0; d<3; d++) {
	  sinks.part[i].acc[d] = force[3*i+d];
	}
      }
  }
    
}

void get_force_and_potential(Bodies & bodies)
{
#ifdef FMM
  if(bodies.size() < 10000)
    {
      direct(bodies);
    }
  else
    {
#if DEBUG
      start("buildTree");
#endif
      Cells cells = buildTree(bodies);

      std::vector<omp_lock_t> p2plocks;
      p2plocks.resize(cells.size());
      std::vector<omp_lock_t> m2llocks;
      m2llocks.resize(cells.size());
      
      //allocate memory for multipoles and local expansions
      std::vector<std::vector<real_t> > Multipoles((int)cells.size(),
						   std::vector<real_t>(NTERM, 0) );

      std::vector<std::vector<real_t>> Locals((int)cells.size(),
					      std::vector<real_t>(NTERM, 0) );
        
      std::vector<std::vector<real_t>> Pns((int)cells.size(),
					   std::vector<real_t>((EXPANSION+1), 0));
                
      for(size_t i=0; i<cells.size(); i++){
	omp_init_lock(&(p2plocks[i]));
        cells[i].p2p_lock = &(p2plocks[i]);

        omp_init_lock(&(m2llocks[i]));
	cells[i].m2l_lock = &(m2llocks[i]);
	
        cells[i].M = &Multipoles[i][0];
        cells[i].L = &Locals[i][0];
        cells[i].Pn= &Pns[i][0];
      }
        
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
	      for(int m =-n; m <= n; m++)
		{
		  int index = n*n + n + m;
		  cells[i].Pn[n] += cells[i].M[index] *
		    cells[i].M[index]*factorial_table[n + m]*factorial_table[n-m];
		}
	      cells[i].Pn[n] = sqrt(cells[i].Pn[n]);
	    }
	}
#if DEBUG
      start("horizontalPass");
#endif

      horizontalPass(cells, cells, true);

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

  direct(bodies);

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
  if(s.last-s.part + 1 != s.n)
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
      if(i >= s.n){
	printf("i:%d s.n:%d dt=%g left=%g right=%g \n", i, s.n, dt, left->timestep, right->timestep);
	ENDRUN("split error 1\n");
      }
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
  slow->n = s.last-left + 1;
  fast->n = (left-s.part);

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
      fast->last = s.part + fast->n-1;
    }

  if(fast->n + slow->n != s.n)
    ENDRUN("split error 2\n");
}

void drift(int clevel, struct sys s, double etime, double dt)
{
#pragma omp parallel for if(s.n > ncrit)
  for(unsigned int i = 0; i < s.n; i++)
    {
      for(int d=0; d<3; d++) {
	COMPSUM(s.part[i].pos[d], s.part[i].pos_e[d], dt*s.part[i].vel[d]);
      }
      s.part[i].postime = etime;
    }
}

void drift_naive(int clevel, struct sys s, double etime)
{
  double dt = etime-s.part[0].postime;
    
#pragma omp parallel for if(s.n > ncrit)
  for(unsigned int i = 0; i < s.n; i++)
    {
      for(int d=0; d<3; d++) {
	COMPSUM(s.part[i].pos[d], s.part[i].pos_e[d], dt*s.part[i].vel[d]);
      }
      s.part[i].postime = etime;
    }
}

void dkd(int clevel, struct sys slow, double stime, double etime, double dt) {
  drift(clevel, slow, stime+dt/2, dt/2);
  if(slow.n>ncrit){
    kick_self(clevel, slow, dt, 1, 0, false);
  }
  else{
    kick_cpu(clevel, slow, slow, dt, 1, 0, false);
  }
  drift(clevel, slow, stime+dt, dt/2);
}

void kdk4(int clevel, struct sys slow, double stime, double etime, double dt) {
  {
      
    if(slow.n>ncrit){
      kick_self(clevel, slow, dt, 1.0/6, 0, false);
    }
    else{
      kick_cpu(clevel, slow, slow, dt, 1.0/6, 0, false);
    }
    
    drift(clevel, slow, stime+dt/2, dt/2);
    
    double b = 0.0, c = 0.0;

    if(slow.n>ncrit){
      kick_self(clevel, slow, dt, b, c, false);
      b=2.0/3.0, c = 1.0/72.0;
      kick_self(clevel, slow, dt, b, c, true);
    }
    else{
      kick_cpu(clevel, slow, slow, dt, b, c, false);
      b=2.0/3.0, c = 1.0/72.0;
      kick_cpu(clevel, slow, slow, dt, b, c, true);
    }
    
    drift(clevel, slow, stime+dt, dt/2);

    if(slow.n>ncrit){
      kick_self(clevel, slow, dt, 1.0/6, 0, false);
    }
    else{
      kick_cpu(clevel, slow, slow, dt, 1.0/6, 0, false);
    }
      
  }
}

/*needs more testing on this integrator...
  void kdk6(int clevel, struct sys slow, double stime, double etime, double dt) {
      
  double b, c;
    
  b = 0.3599508087941436;
  c = 0;
    
  if(slow.n>ncrit){
  kick_self(clevel, slow, dt, b, 0, false);
  }
  else{
  kick_cpu(clevel, slow, slow, dt, b, 0, false);
  }
    
  double a = 1.0798524263824309;
    
  drift(clevel, slow, stime+dt*a, dt*a);
    
  b = -0.1437147273026540;
  c = -0.0139652542242388*0;
    
  if(slow.n>ncrit){
  kick_self(clevel, slow, dt, 0, 0, false);
  kick_self(clevel, slow, dt, b, c, true);
  }
  else{
  kick_cpu(clevel, slow, slow, dt, 0, 0, false);
  kick_cpu(clevel, slow, slow, dt, b, c, true);
  }
    
  a  = -0.5798524263824309;
    
  drift(clevel, slow, etime, dt*a);
    
  b = 0.5675278370170208;
  c = -0.0392470293823456*0;
    
  if(slow.n>ncrit){
  kick_self(clevel, slow, dt, 0, 0, false);
  kick_self(clevel, slow, dt, b, c, true);
  }
  else{
  kick_cpu(clevel, slow, slow, dt, 0, 0, false);
  kick_cpu(clevel, slow, slow, dt, b, c, true);
  }
        
  a  = -0.5798524263824309;
  drift(clevel, slow, etime, dt*a);
    
  b = -0.1437147273026540;
  c = -0.0139652542242388*0;
    
  if(slow.n>ncrit){
  kick_self(clevel, slow, dt, 0, 0, false);
  kick_self(clevel, slow, dt, b, c, true);
  }
  else{
  kick_cpu(clevel, slow, slow, dt, 0, 0, false);
  kick_cpu(clevel, slow, slow, dt, b, c, true);
  }
            
  a = 1.0798524263824309;
      
  drift(clevel, slow, etime, dt*a);
    
  b = 0.3599508087941436;
  c = 0;
    
  if(slow.n>ncrit){
  kick_self(clevel, slow, dt, b, 0, false);
  }
  else{
  kick_cpu(clevel, slow, slow, dt, b, 0, false);
  }
  }
*/

void evolve_split_hold_dkd(int clevel, struct sys total,
                           double stime, double etime,
                           double dt, bool calc_timestep)
{
  struct sys slow = zerosys, fast = zerosys;

  if(calc_timestep) {
    findtimesteps(total);
  }
    
  split(dt, total, &slow, &fast);

  if(fast.n == 0)
    {
      diag->simtime += dt;
#if DEBUG
      printf("t=%g s=%d \n", diag->simtime, total.n);
#endif
    }

  fflush(stdout);

  //hold for fast system
  if(fast.n > 0)
    evolve_split_hold_dkd(clevel+1, fast, stime, stime+dt/2, dt/2, false);

  if(slow.n > 0)
    drift(clevel, slow, stime + dt / 2, dt / 2);
    
  if(slow.n > 0) { //kicksf in between
    if(total.n > ncrit)
      {
	kick_self(clevel, slow, dt, 1, 0, false);

	if(fast.n > 0) {
	  kick_sf(clevel, slow, fast, dt, 1, 0, false);
	  kick_sf(clevel, fast, slow, dt, 1, 0, false);
	}
      }
    else
      {
	kick_cpu(clevel, slow, slow, dt, 1, 0, false);

	if(fast.n > 0) {
	  kick_cpu(clevel, slow, fast, dt, 1, 0, false);
	  kick_cpu(clevel, fast, slow, dt, 1, 0, false);
	}
      }
  }
    
  if(slow.n > 0)
    drift(clevel, slow, etime, dt / 2);

  //hold for fast system
  if(fast.n > 0)
    evolve_split_hold_dkd(clevel+1, fast, stime+dt/2, etime, dt/2, true);
}

void evolve_frost(int clevel, struct sys total,
		  double stime, double etime,
		  double dt, bool calc_timestep)
{
    
  if(total.n == 2) { // change to total.n == 2 for kepler solver
        
    for(int i=0; i<2; i++) {
      evolve_kepler(clevel, total, stime+i*dt/2, stime+(i+1)*dt/2, dt/2);
    }
        
    return;
                    
  }
  // uncomment else out for kepler solver
  else
    {
    
      struct sys slow = zerosys, fast = zerosys;

      if(calc_timestep) {
	findtimesteps(total);
      }
    
      split(dt, total, &slow, &fast);

      if(fast.n == 0)
	{
	  diag->simtime += dt;
#if DEBUG
	  printf("t=%g s=%d \n", diag->simtime, total.n);
#endif
	}

      fflush(stdout);

      if(slow.n > 0 && fast.n > 0) { //kicksf in between
	if(total.n > ncrit)
	  {
	    double b = 1.0/6.0;
	    kick_sf(clevel, slow, fast, dt, b, 0, false);
	    kick_sf(clevel, fast, slow, dt, b, 0, false);
	  }
	else
	  {
	    double b = 1.0/6.0;
	    kick_cpu(clevel, slow, fast, dt, b, 0, false);
	    kick_cpu(clevel, fast, slow, dt, b, 0, false);
	  }
      }
    
    
      //hold for fast system
      if(fast.n > 0)
	evolve_frost(clevel+1, fast, stime, stime+dt/2, dt/2, false);

      //kdk for the slow system
      if(slow.n > 0)
	kdk4(clevel, slow, stime, stime+dt/2, dt/2);
    
      if(slow.n > 0 && fast.n > 0) { //kicksf in between
	if(total.n > ncrit)
	  {
	    double b = 0.0, c = 0.0;
	    kick_sf(clevel, slow, fast, dt, b, c, false);
	    kick_sf(clevel, fast, slow, dt, b, c, false);
	    b = 2.0/3.0, c = 1.0/72.0;
	    kick_sf(clevel, slow, fast, dt, b, c, true);
	    kick_sf(clevel, fast, slow, dt, b, c, true);
	  }
	else
	  {
	    double b = 0.0, c = 0.0;
	    kick_cpu(clevel, slow, fast, dt, b, c, false);
	    kick_cpu(clevel, fast, slow, dt, b, c, false);
	    b = 2.0/3.0, c = 1.0/72.0;
	    kick_cpu(clevel, slow, fast, dt, b, c, true);
	    kick_cpu(clevel, fast, slow, dt, b, c, true);
	  }
      }
    
      //kdk for the slow system
      if(slow.n > 0)
	kdk4(clevel, slow, stime+dt/2, stime+dt, dt/2);
      
      //hold for fast system
      if(fast.n > 0)
	evolve_frost(clevel+1, fast, stime+dt/2, etime, dt/2, true);
    
      if(slow.n > 0 && fast.n > 0) { //kicksf in between
	if(total.n > ncrit)
	  {
	    double b = 1.0/6.0;
	    kick_sf(clevel, slow, fast, dt, b, 0, false);
	    kick_sf(clevel, fast, slow, dt, b, 0, false);
	  }
	else
	  {
	    double b = 1.0/6.0;
	    kick_cpu(clevel, slow, fast, dt, b, 0, false);
	    kick_cpu(clevel, fast, slow, dt, b, 0, false);
	  }
      }
    }
    
}


void do_evolve(struct sys s, double dt)
{
  int clevel;

  if(dt == 0)
    return;
  clevel = 0;

  findtimesteps(s);
    
#if HOLD
  evolve_split_hold_dkd(clevel, s,
			(real_t) diag->simtime,
			(real_t) diag->simtime + dt,
			(real_t) dt, true);
#endif
    
    
#if FROST
  evolve_frost(clevel, s,
	       (real_t) diag->simtime,
	       (real_t) diag->simtime + dt,
	       (real_t) dt, true);
#endif
    
}

double system_kinetic_energy(struct sys s)
{
  unsigned int i;
  long double e = 0.;
  for(i = 0; i < s.n; i++)
    e += s.part[i].mass / G*(s.part[i].vel[0]*s.part[i].vel[0] +
			     s.part[i].vel[1]*s.part[i].vel[1] +
			     s.part[i].vel[2]*s.part[i].vel[2]) / 2;
  return (double) e;
}

double system_potential_energy(struct sys s)
{
  unsigned int i;
  long double e = 0.;
  for(i = 0; i < s.n; i++)
    e += s.part[i].mass / G*s.part[i].pot / 2;
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
	  pos[i] += (long double) s.part[p].mass*s.part[p].pos[i];
	  vel[i] += (long double) s.part[p].mass*s.part[p].vel[i];
	}
      mass += (long double) s.part[p].mass;
    }
  for(int i = 0; i < 3; i++)
    {
      cmpos[i] = pos[i] / mass;
      cmvel[i] = vel[i] / mass;
    }
}

