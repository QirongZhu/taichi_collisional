#include <iostream>
#include <cstring>
#include "exafmm.h"
#include <iostream>
#include <utility>

#include "exafmm.h"
#include "integrator.h"
#include "kernel_fast_lean.h"
#include "timer.h"
#include "traverse_eager.h"
#include "load_gadget.h"
#include "io.h"

using namespace exafmm;

void init_code(void);

double interval;

struct diagnostics global_diag;
struct diagnostics *diag;
struct sys zerosys = { 0, NULL, NULL };

int main(int argc, const char * argv[]) {
  
  P = EXPANSION; //std::stoi(argv[1]);

  //  P = 10;
  initKernel();

  // P2M
  Bodies jbodies(1);
  //  for (int d=0; d<3; d++) jbodies[0].X[d] = 4;
    
  jbodies[0].X[0] = 0;
  jbodies[0].X[1] = 0;
  jbodies[0].X[2] = 0;
  jbodies[0].q    = 1;
  jbodies[0].issource = true;

  jbodies[0].V[0] = 10;
  jbodies[0].V[1] = 10;
  jbodies[0].V[2] = 10;

  Cells cells(4);
  Cell * Cj = &cells[0];
  Cj->X[0] = 0;
  Cj->X[1] = 0;
  Cj->X[2] = 0;

  Cj->R = 1;
  Cj->BODY = &jbodies[0];
  Cj->NBODY = jbodies.size();
  //  Cj->M.resize(NTERM, 0.0);
  //  Cj->Mdot.resize(NTERM, 0.0);

    
  P2M(Cj);

  for(int ind=0; ind<NTERM; ind++)
    printf("M[%d]=%g %g \n", ind, Cj->M[ind], Cj->Mdot[ind]);

  //    printf("GOOD here \n");
  //    fflush(stdout);
    
  // M2M
  Cell * CJ = &cells[1];
  CJ->CHILD = Cj;
  CJ->NCHILD = 1;
  CJ->X[0] = 1;
  CJ->X[1] = 1;
  CJ->X[2] = 2;
  CJ->R = 2;
  CJ->BODY = &jbodies[0];
  CJ->NBODY = jbodies.size();
  //  CJ->M.resize(NTERM, 0.0);
  //  CJ->Mdot.resize(NTERM, 0.0);

  M2M(CJ);
    
  //    printf("GOOD here \n");
  //    fflush(stdout);

  Cell * Ci = &cells[2];
  Ci->X[0] = 120;
  Ci->X[1] = 100;
  Ci->X[2] = 300;
  Ci->R = 2;
    
  //  Ci->L.resize(NTERM, 0.0);
  //  Ci->Ldot.resize(NTERM, 0.0);

  Ci->has_sink=true;

  start("P2M");
  for(int i=0; i<1; i++)
    //    P2M(Cj);
    M2L(Ci, CJ);
  stop("P2M");    

  //  for(int i=0; i<NTERM; i++)
  //  printf("Ci[%d]=%g \n", i, Ci->L[i]);

  Bodies bodies(1);
  
  bodies[0].X[0] = 120;
  bodies[0].X[1] = 100;
  bodies[0].X[2] = 300;
  bodies[0].q = 1;
  bodies[0].p = 0;
    
  bodies[0].issink = true;

  bodies[0].V[0] = 1;
  bodies[0].V[1] = 2;
  bodies[0].V[2] = 3;
    
  for (int d=0; d<3; d++) {
    bodies[0].F[d] = 0;
    bodies[0].J[d] = 0;
  }

  Ci->BODY = &bodies[0];
  Ci->NBODY = bodies.size();
  L2P(Ci);

  // P2P
  Bodies bodies2(1);
  for (size_t b=0; b<bodies2.size(); b++) {
    bodies2[b] = bodies[b];
    bodies2[b].p = 0;
    for (int d=0; d<3; d++) {
      bodies2[b].F[d] = 0;
      bodies2[b].J[d] = 0;
    }
  }
    
  Cj->NBODY = jbodies.size();
  Ci->NBODY = bodies2.size();
  Ci->BODY = &bodies2[0];
    
  P2P(Ci, Cj);

  // Verify results
  real_t pDif = 0, pNrm = 0, FDif = 0, FNrm = 0;
  for (size_t b=0; b<bodies.size(); b++) {
    pDif += (bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p);
    pNrm += bodies[b].p * bodies[b].p;
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]) +
      (bodies[b].F[2] - bodies2[b].F[2]) * (bodies[b].F[2] - bodies2[b].F[2]);
        FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1] +
	  bodies[b].F[2] * bodies[b].F[2];
  }

  printf("Approx force %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e\n",
           bodies[0].p, bodies[0].F[0], bodies[0].F[1], bodies[0].F[2], bodies[0].J[0], bodies[0].J[1], bodies[0].J[2]);
    
    std::cout<<"..........\n"<<std::endl;
    
    //printf("%8.5e %8.5e %8.5e %8.5e\n", bodies2[0].p, bodies2[0].X[0], bodies2[0].X[1], bodies2[0].X[2]);
    
    printf("True force %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e %8.5e \n\n",
           bodies2[0].p, bodies2[0].F[0], bodies2[0].F[1], bodies2[0].F[2],
           bodies2[0].J[0], bodies2[0].J[1], bodies2[0].J[2]);
    
    printf("Diff force %8.5e %8.5e %8.5e %8.5e\n Diff jerk %8.5e %8.5e %8.5e \n",   bodies[0].p-bodies2[0].p,
           bodies[0].F[0]-bodies2[0].F[0],
           bodies[0].F[1]-bodies2[0].F[1],
           bodies[0].F[2]-bodies2[0].F[2],
           bodies[0].J[0]-bodies2[0].J[0],
           bodies[0].J[1]-bodies2[0].J[1],
           bodies[0].J[2]-bodies2[0].J[2]);
}
