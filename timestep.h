#pragma once

#include "exafmm.h"
#include "build_tree.h"

namespace exafmm
{

void timestepCoreTest(Cell * Ci, Cell * Cj)
{
  Body *Bi = Ci->BODY;
  Body *Bj = Cj->BODY;

  int ni = Ci->NBODY;
  int nj = Cj->NBODY;

  for(int i = 0; i < ni; i++)
    {
      real_t timestep = 0, dtimestep = 0;
        
      for(int j = 0; j < nj; j++)
        {
      for(int d = 0; d < 3; d++)
        {
          dX[d] = Bj[j].X[d] - Bi[i].X[d];
          dV[d] = Bj[j].V[d] - Bi[i].V[d];
        }

      real_t R2 = norm(dX);

      if(R2 > 0)
        {
        real_t V2 = norm(dV);
        real_t invR2 = 1/R2;
        real_t invR4 = invR2 * invR2;
        real_t invR6 = invR4 * invR2;
        
        real_t vdotr = (dX[0]*dV[0] + dX[1]*dV[1] + dX[2]*dV[2]);
        real_t msum = (Bi[i].q + Bj[j].q);
            
        real_t Pinv4 = 1.0 / (16.0 * PI4)* (msum*msum) * invR6;
        real_t Qinv4 = 256.0 * (V2 * V2) * invR4;
        timestep += (Qinv4+Pinv4);
            
        real_t dPdt =  3.0 / (32.0 * PI4) *(msum*msum) * invR4 * invR4;
        //real_t dQdt =  256.0 * (V2 + msum * sqrt(invR2)) * sqrt(V2*invR2) * invR6;
        
        real_t vdotr2 = (dX[0]*dV[0]+dX[1]*dV[1]+dX[2]*dV[2]) * invR2;
        real_t invR = sqrt(invR2);
        real_t dQdt = 256.0 * V2 * sqrt(V2) * invR4 * invR * vdotr2 * (1.0 + msum*invR/V2);
            
        dtimestep += (dQdt+dPdt);
                        
        }

#pragma omp atomic
      Bi[i].timestep  += timestep;
#pragma omp atomic
      Bi[i].dtimestep += dtimestep;
    }
}
}

void timestepCore(Cell * Ci, Cell * Cj)
{
  Body *Bi = Ci->BODY;
  Body *Bj = Cj->BODY;

  int ni = Ci->NBODY;
  int nj = Cj->NBODY;

#if SIMD_P2P
  int nii = (ni + NSIMD - 1) & (-NSIMD);
        
#ifndef DOUBLE_P2P
  float Mi[nii] __attribute__ ((aligned(64)));
  float Xi[nii] __attribute__ ((aligned(64)));
  float Yi[nii] __attribute__ ((aligned(64)));
  float Zi[nii] __attribute__ ((aligned(64)));
  float VXi[nii] __attribute__ ((aligned(64)));
  float VYi[nii] __attribute__ ((aligned(64)));
  float VZi[nii] __attribute__ ((aligned(64)));
#else
  double Mi[nii] __attribute__ ((aligned(64)));
  double Xi[nii] __attribute__ ((aligned(64)));
  double Yi[nii] __attribute__ ((aligned(64)));
  double Zi[nii] __attribute__ ((aligned(64)));
  double VXi[nii] __attribute__ ((aligned(64)));
  double VYi[nii] __attribute__ ((aligned(64)));
  double VZi[nii] __attribute__ ((aligned(64)));
#endif

  for(int k = 0; k < ni; k++)
    {
      Mi[k] = Bi[k].q;
      Xi[k] = -Bi[k].X[0];
      Yi[k] = -Bi[k].X[1];
      Zi[k] = -Bi[k].X[2];
      VXi[k] = -Bi[k].V[0];
      VYi[k] = -Bi[k].V[1];
      VZi[k] = -Bi[k].V[2];
    }

#ifndef DOUBLE_P2P
  Vec16f mi, xi, yi, zi, vxi, vyi, vzi, r2, v2, vdotdr2, r, mj;
  Vec16f invR, factor1, fac1, dx, dy, dz, dvx, dvy, dvz;
  Vec16f timestep, tau, dtau;
#else
  Vec8d mi, xi, yi, zi, vxi, vyi, vzi, r2, v2, vdotdr2, r, mj;
  Vec8d invR, factor1, fac1, dx, dy, dz, dvx, dvy, dvz;
  Vec8d timestep, tau, dtau;
#endif
        
  for(int i = 0; i < nii; i = i + NSIMD)
    {
      mi.load(Mi + i);
      xi.load(Xi + i);
      yi.load(Yi + i);
      zi.load(Zi + i);
        
      vxi.load(VXi + i);
      vyi.load(VYi + i);
      vzi.load(VZi + i);

      timestep = 0;

      for(int j = 0; j < nj; j++){
          
      mj = Bj[j].q;
          
      dx = xi + Bj[j].X[0];
      dy = yi + Bj[j].X[1];
      dz = zi + Bj[j].X[2];
          
      dvx = vxi + Bj[j].V[0];
      dvy = vyi + Bj[j].V[1];
      dvz = vzi + Bj[j].V[2];

      r2 = dx * dx + dy * dy + dz * dz;
      v2 = dvx * dvx + dvy * dvy + dvz * dvz + 1e-20;

      mi += mj;
      r2 = select(r2>0, r2, 1e38);

      vdotdr2 = (dx * dvx + dy * dvy + dz * dvz) / r2;

      r = sqrt(r2);
      invR = 1 / r;

      tau = dt_param / M_SQRT2 * r * sqrt(r / mi);
      dtau = 1.5 * tau * vdotdr2;
      dtau = select(dtau<1, 1-dtau/2, 0.5);
      tau = dtau/tau;
      timestep += tau*tau*tau*tau;

      tau = dt_param * r / sqrt(v2);
      dtau = tau * vdotdr2 * (1 + mi / (v2 * r));
      dtau = select(dtau<1, 1-dtau/2, 0.5);
      tau = dtau/tau;

      timestep += tau*tau*tau*tau;
      }
        
      timestep.store(Xi+i);
    }
        
    for(int i = 0; i < ni; i++) {
#pragma omp atomic
    Bi[i].timestep += (real_t) Xi[i];
    }
#else
  for(int i = 0; i < ni; i++)
    {
      real_t timestep = 0;

      for(int j = 0; j < nj; j++)
        {
      for(int d = 0; d < 3; d++)
        {
          dX[d] = Bj[j].X[d] - Bi[i].X[d];
          dV[d] = Bj[j].V[d] - Bi[i].V[d];
        }

      real_t R2 = norm(dX);

      if(R2 > 0) {
          real_t R = sqrt(R2);
          
          real_t vdotdr2;
          
          real_t v2 = norm(dV) + 1e-20;

          vdotdr2 = (dX[0] * dV[0] + dX[1] * dV[1] + dX[2] * dV[2]) / R2;

          real_t tau = dt_param / M_SQRT2 * sqrt(R * R2 / (Bi[i].q + Bj[j].q));
          real_t dtau = 3 * tau * vdotdr2 / 2;
          
          if(dtau > 1)
            dtau = 1;
          tau = (1 - dtau / 2)/tau;
        
          timestep += tau*tau*tau*tau;

          if(v2 > 0) {
              tau = dt_param * R / sqrt(v2);
              dtau = tau * vdotdr2 * (1 + (Bi[i].q + Bj[j].q) / (v2 * R));
              if(dtau > 1)
                  dtau = 1;
              tau = (1 - dtau / 2)/tau;
              timestep += tau*tau*tau*tau;
          }
        }
            
#pragma omp atomic
      Bi[i].timestep += timestep;
      
    }
    
#endif
}

  //! Recursive call to dual tree traversal for horizontal pass
void timestepPass(Cell * Ci, Cell * Cj)
  {
    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];	// Distance vector from source to target

    real_t R2 = norm(dX) * theta * theta;	// Scalar distance squared

    if(R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R)) { // If distance is far enough
        if((Ci->NBODY <= 4 && Cj->NCHILD == 0) ||
           (Cj->NBODY <= 4 && Ci->NCHILD == 0)) {
                timestepCore(Ci, Cj);
        }
        else{
            return;
        }
    }
    else if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {				// Else if both cells are leafs
          timestepCore(Ci, Cj);
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {				// If Cj is leaf or Ci is larger
	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {			// Loop over Ci's children
#pragma omp task untied if(ci->NBODY > 100)	//   Start OpenMP task if large enough task
          timestepPass(ci, Cj);	//   Recursive call to target child cells
	  }			//  End loop over Ci's children
      }
    else
      {				// Else if Ci is leaf or Cj is larger
          for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++){// Loop over Cj's children
          timestepPass(Ci, cj);	//   Recursive call to source child cells
          }			//  End loop over Cj's children
      }				// End if for leafs and Ci Cj size
  }

  //! Horizontal pass interface
  void timestepPass(Cells & icells)
  {
#pragma omp parallel
#pragma omp single nowait
      timestepPass(&icells[0], &icells[0]);
  }

    
void findtimesteps(struct sys system) {

    if(system.n > ncrit/2) {
    
    bodies.clear();
    bodies.resize(system.n);

    for(unsigned int i = 0; i < system.n; i++)
      {
        bodies[i].X = system.part[i].pos;
        bodies[i].V = system.part[i].vel;
        bodies[i].index = i;
        bodies[i].q = system.part[i].mass;
        bodies[i].timestep = 0;
        bodies[i].dtimestep = 0;
      }
    
    Cells cells = buildTree(bodies);

    timestepPass(cells);
    
#pragma omp parallel for if(bodies.size() > ncrit)
    for(size_t b = 0; b < bodies.size(); b++) {
        unsigned int i = bodies[b].index;
        real_t tau    = 1/sqrt(sqrt(bodies[b].timestep));
        system.part[i].timestep = tau;
      }
    }
    else{
        
        real_t min_step = HUGE;
        
        for(unsigned int i = 0; i < system.n; i++) {
            
        real_t timestep = 0;
              
        for(unsigned int j = 0; j < system.n; j++) {
            
             if(system.part[i].id == system.part[j].id)
               continue;

                for(int d=0; d<3; d++){
                    dX[d] = system.part[j].pos[d] - system.part[i].pos[d];
                    dV[d] = system.part[j].vel[d] - system.part[i].vel[d];
                }
                
                real_t R2 = norm(dX);

                real_t R = sqrt(R2);
                real_t vdotdr2;
                real_t v2 = norm(dV) + 1e-20;

                vdotdr2 = (dX[0] * dV[0] + dX[1] * dV[1] + dX[2] * dV[2]) / R2;

                real_t tau = dt_param / M_SQRT2 * sqrt(R * R2 / (system.part[i].mass + system.part[j].mass));
                real_t dtau = 3 * tau * vdotdr2 / 2;
                        
                if(dtau > 1) dtau = 1;
                tau = (1 - dtau / 2)/tau;
                timestep += tau*tau*tau*tau;

                tau = dt_param * R / sqrt(v2);
                dtau = tau * vdotdr2 * (1 + (system.part[i].mass + system.part[j].mass) / (v2 * R));

                if(dtau > 1) dtau = 1;
                tau = (1 - dtau / 2)/tau;
                timestep += tau*tau*tau*tau;
              }
            
            system.part[i].timestep = 1.0/sqrt(sqrt(timestep));
                        
            //if(min_step > system.part[i].timestep)
            //    min_step = system.part[i].timestep;
          }
                
        //        for(unsigned int i = 0; i < system.n; i++) {
        //            system.part[i].timestep = min_step;
        //        }
        
    }
}


}
