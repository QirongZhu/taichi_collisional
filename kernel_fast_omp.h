#pragma once

#include "exafmm.h"
#include <iostream>
#include <math.h>        /*isnan, sqrt, atan2 */
#include <stdio.h>        /*printf */

#include "vectorclass/vectorclass.h"
#ifndef DOUBLE_P2P
#define NSIMD 16
#else
#define NSIMD 8
#endif

#ifdef MINIBALL
#include "miniball.hpp"
#endif

#include "kernel_rotate.h"

namespace exafmm
{
  real_t oddOrEventable[EXPANSION+1];

  real_t factorial_table[2*EXPANSION+1];
  real_t inverse_factorial_table[2*EXPANSION+1];

  real_t factorial_coef[EXPANSION+1][EXPANSION+1];
  real_t factorial_coef_inv[EXPANSION+1][EXPANSION+1];

  real_t factorial_coef_oned[(EXPANSION+1)*(EXPANSION+1)];
  real_t factorial_coef_inv_oned[(EXPANSION+1)*(EXPANSION+1)];

  real_t combinator_coef[EXPANSION+1][EXPANSION+1];

  const complex_t R(1., 0.);    //!<Real unit
  const complex_t I(0., 1.);    //!<Imaginary unit

  real_t norm(real_t*X)
  {
    return X[0]*X[0]+X[1]*X[1]+X[2]*X[2];
  }


  int oddOrEven(int n)
  {
    return (((n) & 1) == 1) ? -1 : 1;
  }

  void initKernel()
  {
    NTERM = (P+1)*(P+1);
    factorial_table[0] = 1.0;
    inverse_factorial_table[0] = 1.0;
    for(int i = 1; i < EXPANSION*2; i++)
      {
	factorial_table[i] = (real_t) i *factorial_table[i - 1];
	inverse_factorial_table[i] = 1 / factorial_table[i];
      }
    for(int i = 0; i <= P; i++)
      {
	for(int j = 0; j <= P; j++)
	  {
	    factorial_coef[i][j] = factorial_table[i+j]*factorial_table[i - j];
	    factorial_coef_inv[i][j] = 1.0 / factorial_coef[i][j];
	    if(j <= i)
	      combinator_coef[i][j] =
		factorial_table[i]*inverse_factorial_table[j]*inverse_factorial_table[i - j];
	  }

	oddOrEventable[i] = oddOrEven(i);

      }

    for(int n = 0; n <= P; n++)
      {
	int index_start = n*n+n;
	factorial_coef_oned[index_start] = factorial_coef[n][0];
	factorial_coef_inv_oned[index_start] = factorial_coef_inv[n][0];
	for(int m = 1; m <= n; m++)
	  {
	    factorial_coef_oned[index_start+m] = factorial_coef[n][m];
	    factorial_coef_oned[index_start - m] = factorial_coef[n][m];
	    factorial_coef_inv_oned[index_start+m] = factorial_coef_inv[n][m];
	    factorial_coef_inv_oned[index_start - m] = factorial_coef_inv[n][m];
	  }
      }
  }

  void real_2_complex(real_t*real_arr, complex_t*complex_arr, int order)
  {
    for(int n = 0; n <= order; n++)
      {
	int index_start = n*n+n;
	complex_arr[index_start] = complex_t(real_arr[index_start], 0);
	real_t oddevenfac = -1;
	for(int m = 1; m <= n; m++)
	  {
	    complex_arr[index_start - m] = complex_t(real_arr[index_start+m]*oddevenfac,
						     -real_arr[index_start - m]*oddevenfac);
	    complex_arr[index_start+m] = complex_t(real_arr[index_start+m], real_arr[index_start - m]);
	    oddevenfac *= -1;
	  }
      }
  }

  void complex_2_real(complex_t*complex_arr, real_t*real_arr, int order)
  {
    for(int n = 0; n <= order; n++)
      {
	int index_start = n*n+n;
	real_arr[index_start] = std::real(complex_arr[index_start]);
	for(int m = 1; m <= n; m++)
	  {
	    real_arr[index_start - m] = std::imag(complex_arr[index_start+m]);
	    real_arr[index_start+m] = std::real(complex_arr[index_start+m]);
	  }
      }
  }

  void make_Tnm(real_t*dX, complex_t*Tnm, int Order)
  {
    real_t r2 = norm(dX);
    real_t invr2 = 1 / r2;
    Tnm[0] = complex_t(sqrt(invr2), 0);

    for(int n = 1; n <= Order; n++)
      {
	int n2 = n*n;
	int ind = n2+n+n;

	Tnm[ind] = Tnm[n2 - 1]*complex_t(dX[0]*(2*n - 1)*invr2, dX[1]*(2*n - 1)*invr2);

	if(n > 0)
	  {
	    int m, indexstart = n2 - n, indexnew = n2+n;

	    if(n > 1)
	      {
		for(m = 0; m < n - 1; m++)
		  {
		    Tnm[indexnew+m] =
		      invr2*(2*n - 1.)*dX[2]*Tnm[indexstart+m] -
		      invr2*((n - 1)*(n - 1)*(n > 1) - m*m)*Tnm[n2 - 3*n+2+m];
		  }
	      }

	    m = n - 1;
	    Tnm[indexnew+m] = invr2*(2*n - 1.)*dX[2]*Tnm[indexstart+m];
	  }

	real_t oddoreven = -1;

	for(int m = 1; m <= n; m++)
	  {
	    Tnm[n2+n - m] = std::conj(Tnm[n2+n+m]*oddoreven);
	    oddoreven *= -1;
	  }
      }
  }

  void make_Gnm(real_t*dX, complex_t*Gnm, int Order)
  {
    real_t r2 = norm(dX);

    Gnm[0] = complex_t(1, 0);

    for(int n = 1; n <= Order; n++)
      {
	int n2 = n*n;
	int ind = n2+n+n;

	Gnm[ind] = Gnm[n2 - 1]*complex_t(dX[0] / (2*n), dX[1] / (2*n));

	int m, m2, indexstart = n2 - n, indexnew = n2+n;

	for(m = 0; m < n - 1; m++)
	  {
	    m2 = n2 - m*m;
	    Gnm[indexnew+m] =
	      (2*n - 1.) / m2*dX[2]*Gnm[indexstart+m] - r2 / m2*Gnm[n2 - 3*n+2+m];
	  }
        
	Gnm[indexnew+n - 1] = dX[2]*Gnm[indexstart+n - 1];

	real_t oddoreven = -1;
	for(int m = 1; m <= n; m++)
	  {
	    Gnm[indexnew - m] = std::conj(Gnm[indexnew+m]*oddoreven);
	    oddoreven *= -1;
	  }
      }
  }

void make_Gnm_conj(real_t*dX, complex_t*Gnm, int Order)
{
  real_t r2 = norm(dX);

  Gnm[0] = complex_t(1, 0);

  for(int n = 1; n <= Order; n++)
    {
  int n2 = n*n;
  int ind = n2+n+n;

  Gnm[ind] = Gnm[n2 - 1]*complex_t(dX[0] / (2*n), -dX[1] / (2*n));

  int m, m2, indexstart = n2 - n, indexnew = n2+n;

  for(m = 0; m < n - 1; m++)
    {
      m2 = n2 - m*m;
      Gnm[indexnew+m] =
        (2*n - 1.) / m2*dX[2]*Gnm[indexstart+m] - r2 / m2*Gnm[n2 - 3*n+2+m];
    }
      
  Gnm[indexnew+n - 1] = dX[2]*Gnm[indexstart+n - 1];

  real_t oddoreven = -1;
  for(int m = 1; m <= n; m++)
    {
      Gnm[indexnew - m] = std::conj(Gnm[indexnew+m]*oddoreven);
      oddoreven *= -1;
    }
    }
}


  void P2P(Cell*Ci, Cell*Cj)
  {
      
    if(Ci->has_sink)
      {
	Body *Bi = Ci->BODY;
	Body *Bj = Cj->BODY;

	int ni = Ci->NBODY;
	int nj = Cj->NBODY;

#if SIMD_P2P
	if(ni > NSIMD)
	  {
	    int nii = (ni+NSIMD - 1) & (-NSIMD);

#ifndef DOUBLE_P2P
	    float Xi[nii],Yi[nii],Zi[nii],Mi[nii];
#else
	    double Xi[nii],Yi[nii],Zi[nii],Mi[nii];
#endif
          
	    for(int k = 0; k < ni; k++)
	      {
		Xi[k]  = -Bi[k].X[0];
		Yi[k]  = -Bi[k].X[1];
		Zi[k]  = -Bi[k].X[2];
	      }
          
#ifndef DOUBLE_P2P
	    Vec16f xi, yi, zi, r, mj, mi;
	    Vec16f invR, dx, dy, dz, r2;
	    Vec16f ax, ay, az, pot;
#else
	    Vec8d xi, yi, zi, r, mj, mi;
	    Vec8d invR, dx, dy, dz, r2;
	    Vec8d ax, ay, az, pot;
#endif
	    for(int i = 0; i < nii; i = i+NSIMD)
	      {
		xi.load(Xi+i);
		yi.load(Yi+i);
		zi.load(Zi+i);
              
		ax = 0, ay = 0, az = 0, pot = 0;

		for(int j = 0; j < nj; j++)
		  {
            if(!Bj[j].issource)
                continue;

		    mj = Bj[j].q;
		    dx = Bj[j].X[0]+xi;
		    dy = Bj[j].X[1]+yi;
		    dz = Bj[j].X[2]+zi;
		    r2 = dx*dx+dy*dy+dz*dz;

		    r    = sqrt(r2);
		    invR = 1 / r;
            invR = select(r2 > 0, invR, 0);
		    mj  *= invR;
              
		    pot += mj;
              
		    mj  = mj*(invR*invR);
		    ax += dx*mj;
		    ay += dy*mj;
		    az += dz*mj;
		  }

		pot.store(Mi+i);
		ax.store(Xi+i);
		ay.store(Yi+i);
		az.store(Zi+i);
	      }
          
	    for(int i = 0; i < ni; i++) {
        if(Bi[i].issink)
		{
		  #pragma omp atomic
		  Bi[i].p    += (real_t) Mi[i];
		  #pragma omp atomic
		  Bi[i].F[0] += (real_t) Xi[i];
		  #pragma omp atomic
		  Bi[i].F[1] += (real_t) Yi[i];
		  #pragma omp atomic
		  Bi[i].F[2] += (real_t) Zi[i];
                }
            }
          
	  }
	else
	  {
	    // do not use vectorized version if n1*n2 too small
	    for(int i = 0; i < ni; i++)
	      {
		if(!Bi[i].issink)
		  continue;
        
		real_t acc[3] = {0.0,0.0,0.0};
		real_t dx[3] = {0.0,0.0,0.0};

		real_t pot = 0;

		for(int j = 0; j < nj; j++)
		  {
              if(!Bj[j].issource)
                    continue;

		    for(int d=0; d<3; d++)
		      dx[d] = Bj[j].X[d] - Bi[i].X[d];

		    real_t R2 = (dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

		    if(R2 > 0)
            {
		      real_t invR2 = 1 / R2;
		      real_t invR = Bj[j].q*sqrt(invR2);
		      pot += invR;
		      for(int d=0; d<3; d++){
			dx[d] *= invR2*invR;
			acc[d] += dx[d];
		      }
		    }
		  }
		#pragma omp atomic
		Bi[i].p += (real_t) pot;
		#pragma omp atomic
		Bi[i].F[0] += (real_t) acc[0];
		#pragma omp atomic
		Bi[i].F[1] += (real_t) acc[1];
		#pragma omp atomic
		Bi[i].F[2] += (real_t) acc[2];
	      }
	  }
#else
	// do not use vectorized version if n1*n2 too small
	for(int i = 0; i < ni; i++)
	  {
	    real_t acc[3] = {0.0,0.0,0.0};
	    real_t dx[3] = {0.0,0.0,0.0};
          
	    real_t pot = 0;

	    for(int j = 0; j < nj; j++)
	      {
            if(!Bj[j].issource)
                continue;

		for(int d=0; d<3; d++)
		  dx[d] = Bi[i].X[d] - Bj[j].X[d];

		real_t R2 = norm(dX);

		if(R2 > 0)
		  {
		    real_t R = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
		    real_t invR2 = 1.0 / R2;
		    real_t invR = Bj[j].q*sqrt(invR2) ;
		    pot += invR;
              
		    for(int d=0; d<3; d++) {
		      dx[d] *= invR2*invR;
		      acc[d] += dx[d];
		    }
		  }
	      }

	    if(Bi[i].issink) {
	      #pragma omp atomic
	      Bi[i].p += pot;
	      #pragma omp atomic
	      Bi[i].F[0] -= acc[0];
	      #pragma omp atomic
	      Bi[i].F[1] -= acc[1];
	      #pragma omp atomic
	      Bi[i].F[2] -= acc[2];
	    }
	  }
#endif

      }
  }

  void P2M(Cell*C)
  {
    real_t min_acc = 1e30;

#ifndef MINIBALL
    if(C->NBODY > 0)
      {
	real_t comx = 0, comy = 0, comz = 0, totalq = 0;

	for(Body*B = C->BODY; B != C->BODY+C->NBODY; B++)
	  {
	    comx += B->X[0]*B->q;
	    comy += B->X[1]*B->q;
	    comz += B->X[2]*B->q;
	    totalq += B->q;
	  }

	comx /= totalq;
	comy /= totalq;
	comz /= totalq;

	C->X[0] = comx;
	C->X[1] = comy;
	C->X[2] = comz;
	C->cell_mass = totalq;
      }
    else
      {
	C->R = 1e-10;
      }
#else
    if(C->NBODY > 0)
      {
	std::list < std::vector < real_t > >lp;

        for(Body*B = C->BODY; B != C->BODY+C->NBODY; B++)
	  {
	    std::vector < real_t > p = {B->X[0],  B->X[1],  B->X[2]};
	    lp.push_back(p);
	  }
          
	// define the types of iterators through the points and their coordinates
	typedef std::list < std::vector < real_t > >::const_iterator PointIterator;
	typedef std::vector < real_t >::const_iterator CoordIterator;
	typedef Miniball::Miniball < Miniball::CoordAccessor < PointIterator, CoordIterator > >MB;
	MB mb(3, lp.begin(), lp.end());
	real_t rminball = sqrt(mb.squared_radius());
	const real_t *center = mb.center();
          
	//use rminball as R
	C->X[0] = center[0];
	C->X[1] = center[1];
	C->X[2] = center[2];
	C->R = rminball+1e-10;

      }
    else
      {
	C->R = 1e-10;
      }
#endif

    complex_t c_multipole[NTERM];

#ifndef MINIBALL
    max_r2  = 1e-6*C->R*C->R;
#endif

    for(Body*B = C->BODY; B != C->BODY+C->NBODY; B++)
      {
	for(int d = 0; d < 3; d++)
	  {
	    dX[d] = B->X[d] - C->X[d];
	  }

#ifndef MINIBALL
	max_r2 = std::max(max_r2, norm(dX));
#endif

	min_acc = std::min(B->acc_old, min_acc);

	if(B->issink)
	  C->has_sink = true;

	complex_t Gnm[NTERM];

	make_Gnm(&dX[0], &Gnm[0], P);

	for(int indice = 0; indice < NTERM; indice++)
	  c_multipole[indice] += B->q*Gnm[indice];
      }

    real_t r_multipole[NTERM];
    complex_2_real(&c_multipole[0], &r_multipole[0], P);

    for(int indice = 0; indice < NTERM; indice++)
      C->M[indice] += r_multipole[indice];

    C->min_acc = min_acc;

#ifndef MINIBALL
    C->R       = sqrt(max_r2);
#endif

  }

  void M2M(Cell*Ci)
  {
    real_t min_acc = 1e30;

#ifndef MINIBALL
    if(Ci->NCHILD > 0)
      {
	real_t comx = 0, comy = 0, comz = 0, totalq = 0;

	for(Cell*Cj = Ci->CHILD; Cj != Ci->CHILD+Ci->NCHILD; Cj++)
	  {
	    comx += Cj->X[0]*Cj->cell_mass;
	    comy += Cj->X[1]*Cj->cell_mass;
	    comz += Cj->X[2]*Cj->cell_mass;
	    totalq += Cj->cell_mass;
	    min_acc = std::min(Cj->min_acc, min_acc);

	    if(Cj->has_sink)
	      Ci->has_sink = true;
	  }

	comx /= totalq;
	comy /= totalq;
	comz /= totalq;

	real_t max_r = 1e-10;

	for(Cell*Cj = Ci->CHILD; Cj != Ci->CHILD+Ci->NCHILD; Cj++)
	  {
	    for(int d = 0; d < 3; d++)
	      dX[d] = Cj->X[d] - Ci->X[d];
	    max_r = std::max(max_r, sqrt(norm(dX))+Cj->R);
	  }

	Ci->R = max_r;
	Ci->min_acc = min_acc;
	Ci->X[0] = comx;
	Ci->X[1] = comy;
	Ci->X[2] = comz;
	Ci->cell_mass = totalq;
      }
    else
      {
	Ci->R = 1e-10;
	Ci->min_acc = min_acc;
      }
#else
    if(Ci->NCHILD > 0)
      {

	std::vector<Ball> S;

	for(Cell*ci = Ci->CHILD; ci != Ci->CHILD+Ci->NCHILD; ci++)
	  {
	    min_acc = std::min(ci->min_acc, min_acc);

	    if(ci->has_sink)
	      Ci->has_sink = true;

	    if(ci->NCHILD == 0)
	      {            //leaves, just use bodies to get R_max
		    double a[3] = {ci->X[0], ci->X[1], ci->X[2]};
		    S.push_back(Ball(a, ci->R));
	      }
	    else
	      { // use granddaughers information to get R_max
              for(Cell*cii = ci->CHILD; cii != ci->CHILD+ci->NCHILD; cii++)
              {
                  double a[3] = {cii->X[0], cii->X[1], cii->X[2]};
                  S.push_back(Ball(a, cii->R));
              }
	      }
	  }
          
        Minsphere mb(S.begin(), S.end());

        real_t rminball = mb.radius();

	const real_t *center = mb.center_cartesian_begin();
        Ci->X[0] = center[0];
        Ci->X[1] = center[1];
        Ci->X[2] = center[2];
        Ci->R    = rminball ;
          
	Ci->min_acc = min_acc;
      }
    else
      {
	Ci->R = 1e-10;
	Ci->min_acc = min_acc;
      }
#endif

    complex_t c_multipole[NTERM];
    for(int i = 0; i < NTERM; i++)
      {
	c_multipole[i] = 0;
      }

    for(Cell*Cj = Ci->CHILD; Cj != Ci->CHILD+Ci->NCHILD; Cj++)
      {
	for(int d = 0; d < 3; d++)
	  {
	    dX[d] = Cj->X[d] - Ci->X[d];
	  }

	real_t r_multipole_Cj[NTERM];
	complex_t c_multipole_Cj[NTERM];

	for(int i = 0; i < NTERM; i++)
	  {
	    r_multipole_Cj[i] = Cj->M[i];
	    c_multipole_Cj[i] = 0;
	  }

	real_2_complex(r_multipole_Cj, c_multipole_Cj, P);

	complex_t Gnm[NTERM];
	make_Gnm(&dX[0], &Gnm[0], P);

	for(int n = 0; n <= P; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= n; k++)
		  for(int l = std::max(-k, m - n+k); l <= std::min(k, m+n - k); l++)
		    c_multipole[indice] +=
		      c_multipole_Cj[(n - k)*(n - k)+(n - k)+(m - l)]*Gnm[k*k+k+l];
	      }
	  }
      }

    real_t r_multipole[NTERM];
    complex_2_real(c_multipole, r_multipole, P);

    for(int indice = 0; indice < NTERM; indice++)
      Ci->M[indice] += r_multipole[indice];
  }

  void L2L(Cell*Ci)
  {
    real_t r_local_Ci[NTERM];
    complex_t c_local_Ci[NTERM];

    for(int indice = 0; indice < NTERM; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, P);

    for(Cell*Cj = Ci->CHILD; Cj != Ci->CHILD+Ci->NCHILD; Cj++)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - Cj->X[d];

	complex_t c_local[NTERM];

	for(int indice = 0; indice < NTERM; indice++)
	  c_local[indice] = 0.0;

	complex_t Gnm[NTERM];

	make_Gnm_conj(&dX[0], &Gnm[0], P);

	for(int n = 0; n <= P; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= P - n; k++)
		  for(int l = -k; l <= k; l++)
		    c_local[indice] += c_local_Ci[(n+k)*(n+k)+(n+k)+(m+l)]*Gnm[k*k+k+l];
	      }
	  }

	real_t r_local[NTERM];

	complex_2_real(c_local, r_local, P);

	for(int indice = 0; indice < NTERM; indice++)
	  Cj->L[indice] += r_local[indice];

      }
  }

  void L2P(Cell*Ci)
  {
    real_t r_local_Ci[NTERM];
    complex_t c_local_Ci[NTERM];

    for(int indice = 0; indice < NTERM; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, P);

    for(Body*B = Ci->BODY; B != Ci->BODY+Ci->NBODY; B++)
      {
	complex_t Phi[4];

	if(!B->issink)
	  continue;

	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - B->X[d];

	complex_t Gnm[NTERM];
          
	make_Gnm_conj(&dX[0], &Gnm[0], P);

	for(int n = 0; n <= 1; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;

		for(int k = 0; k <= P - n; k++)
		  for(int l = -k; l <= k; l++)
		    Phi[indice] += c_local_Ci[(n+k)*(n+k)+(n+k)+(m+l)]*Gnm[k*k+k+l];
	      }
	  }

	B->p    += std::real(Phi[0]);
	B->F[0] -= std::real(Phi[3]);
	B->F[1] -= std::imag(Phi[3]);
	B->F[2] -= std::real(Phi[2]);
      }
  }

  void P2L(Cell*Ci, Cell*Cj)
  {
  
    for (Body*B = Cj->BODY; B != Cj->BODY+Cj->NBODY; B++) {
      for(int d = 0; d < 3; d++)
        dX[d] = Ci->X[d] - B->X[d];

      real_t r_local[NTERM];
      complex_t c_local[NTERM];
      complex_t Tnm[NTERM];

      make_Tnm(&dX[0], &Tnm[0], P);
      for(int indice = 0; indice < NTERM; indice++)
	c_local[indice] = Tnm[indice]*B->q;

      complex_2_real(c_local, r_local, P);

      for(int indice = 0; indice < NTERM; indice++){
	#pragma omp atomic
	Ci->L[indice] += r_local[indice];
      }
    }
  
  }


  void M2P(Cell*Ci, Cell*Cj)
  {
  
    real_t r_multipole_Cj[NTERM];
    complex_t c_Cj[NTERM];

    for(int index = 0; index < NTERM; index++)
      {
	r_multipole_Cj[index] = Cj->M[index];
      }

    for(int n = 0; n <= P; n++)
      {
	int index_start = n*n+n;
	for(int m = 1; m <= n; m++)
	  r_multipole_Cj[index_start - m] *= -1;
      }
  
    real_2_complex(r_multipole_Cj, c_Cj, P);

    for(Body*B = Ci->BODY; B != Ci->BODY+Ci->NBODY; B++)
      {
	complex_t Phi[4];

        for(int d = 0; d < 3; d++)
	  dX[d] = B->X[d] - Cj->X[d];

	std::vector < complex_t > Tnm(NTERM, 0);

	make_Tnm(&dX[0], &Tnm[0], P);

	for(int n = 0; n <= 1; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= P - n; k++)
		  {
		    for(int l = -k; l <= k; l++)
		      {
			int indice2 = k*k+k+l;
			Phi[indice] += c_Cj[indice2]*Tnm[(n+k)*(n+k)+(n+k)+(m+l)];
		      }
		  }
	      }
	  }
	#pragma omp atomic
	B->p    += std::real(Phi[0]);
	#pragma omp atomic
	B->F[0] -= std::real(Phi[3]);
	#pragma omp atomic
	B->F[1] -= std::imag(Phi[3]);
	#pragma omp atomic
	B->F[2] -= std::real(Phi[2]);
      }
  }

  void M2L_rotate(Cell*Ci, Cell*Cj)
  {
      
    if(Ci->NBODY*Cj->NBODY <= 5000){
        P2P(Ci, Cj);
        return;
    }
      
    if(Ci->NBODY <= 8){
      M2P(Ci, Cj);
      return;
    }
      
    if(Cj->NBODY <= 8)
      {
	P2L(Ci, Cj);
	return;
      }

    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];

    real_t r2_xy = dX[0]*dX[0]+dX[1]*dX[1];
    real_t r_xy  = sqrt(r2_xy);
    real_t r2    = r2_xy+dX[2]*dX[2];
    real_t r     = sqrt(r2);
    real_t r_inv = 1.0 / r;

    real_t exp_a_z_re[P+1], exp_a_z_im[P+1];
    real_t exp_a_x_re[P+1], exp_a_x_im[P+1];

    exp_a_z_re[0] = 1;
    exp_a_z_im[0] = 0;
    exp_a_x_re[0] = 1;
    exp_a_x_im[0] = 0;
    exp_a_z_re[1] = 1;
    exp_a_z_im[1] = 0;
    exp_a_x_re[1] = dX[2]*r_inv;
    exp_a_x_im[1] = -r_xy*r_inv;

    if(r_xy > 0)
      {
	r_xy = 1.0 / r_xy;
	exp_a_z_re[1] = dX[1]*r_xy;
	exp_a_z_im[1] = dX[0]*r_xy;
      }

    for(int m = 1; m < P; m++)
      {
	exp_a_z_re[m+1] = exp_a_z_re[m]*exp_a_z_re[1] - exp_a_z_im[m]*exp_a_z_im[1];
	exp_a_z_im[m+1] = exp_a_z_re[m]*exp_a_z_im[1] + exp_a_z_im[m]*exp_a_z_re[1];
	exp_a_x_re[m+1] = exp_a_x_re[m]*exp_a_x_re[1] - exp_a_x_im[m]*exp_a_x_im[1];
	exp_a_x_im[m+1] = exp_a_x_re[m]*exp_a_x_im[1] + exp_a_x_im[m]*exp_a_x_re[1];
      }

    real_t r_multipole_Cj[NTERM], r_local_ci[NTERM];

    //Step 0:       copy multipoles, and make multipoles into homogenius polynomials
    for(int indice = 0; indice < NTERM; indice++)
      {
	r_local_ci[indice] = 0;
	r_multipole_Cj[indice] = Cj->M[indice]*factorial_coef_oned[indice];
      }

    //Step 1:       first rotate multipoles around the original z - axis with alpha_z
    for(int n = 1; n <= P; n++)
      {
	int index_start = n*n+n;
	real_t Re_ei, Im_ei, Im_r, Re_r;
	real_t oddoreven = -1;
	for(int m = 1; m <= n; m++)
	  {
	    Re_ei = exp_a_z_re[m];
	    Im_ei = exp_a_z_im[m];
	    Im_r = r_multipole_Cj[index_start - m];
	    Re_r = r_multipole_Cj[index_start+m];
	    r_multipole_Cj[index_start - m] = Re_ei*Im_r+Im_ei*Re_r;
	    r_multipole_Cj[index_start+m] = Re_ei*Re_r - Im_ei*Im_r;
	  }
      }

    //Step 2:       swap x and z
    swap_x_z(r_multipole_Cj);

    //step 3:       rotate around(new) z - axis by alpha_x
    for(int n = 1; n <= P; n++)
      {
	int index_start = n*n+n;
	real_t Re_ei, Im_ei, Im_r, Re_r;
	real_t oddoreven = -1;
	for(int m = 1; m <= n; m++)
	  {
	    Re_ei = exp_a_x_re[m];
	    Im_ei = exp_a_x_im[m];
	    Im_r = r_multipole_Cj[index_start - m];
	    Re_r = r_multipole_Cj[index_start+m];
	    r_multipole_Cj[index_start - m] = Re_ei*Im_r+Im_ei*Re_r;
	    r_multipole_Cj[index_start+m] = Re_ei*Re_r - Im_ei*Im_r;
	  }
      }

    //Step 4:       swap x and z
    swap_x_z(r_multipole_Cj);

    //Lastly multiply multipoles by proper normalizations
    for(int n = 0; n <= NTERM; n++)
      {
	r_multipole_Cj[n] *= factorial_coef_inv_oned[n];
      }

    real_t powers_of_r, powers_of_r_n, powers_of_r_m;

    powers_of_r_n = r_inv;

    for(int n = 0; n <= P; n++)
      {
	int index_start = n*(n+1);
	powers_of_r = powers_of_r_n;

	for(int k = 0; k <= P - n; k++)
	  {
	    real_t fac = factorial_table[n+k]*powers_of_r;
	    r_local_ci[index_start] += r_multipole_Cj[k*(k+1)]*fac;
	    powers_of_r *= r_inv;
	  }

	powers_of_r_m = r_inv;

	real_t oddevenfac = -1, fac;
	int index_start_k;

	for(int m = 1; m <= std::min(n, P - n); m++)
	  {
	    powers_of_r = powers_of_r_n*powers_of_r_m;
	    for(int k = m; k <= P - n; k++)
	      {
                index_start_k = k*(k+1);
		fac = factorial_table[n+k]*powers_of_r*oddevenfac;
		r_local_ci[index_start - m] += r_multipole_Cj[index_start_k - m]*fac;
		r_local_ci[index_start+m] += r_multipole_Cj[index_start_k+m]*fac;
		powers_of_r *= r_inv;
	      }
	    oddevenfac *= -1;
	    powers_of_r_m *= r_inv;
	  }
	powers_of_r_n *= r_inv;
      }

    //Step 5:       swap x and z
    swap_x_z(r_local_ci);

    //Step 6:       Rotate local expansion around z - axis by(-alpha_x)
    for(int n = 1; n <= P; n++)
      {
	int index_start = n*n+n;
	real_t Re_ei, Im_ei, Im_c, Re_c;
	real_t oddoreven = -1;
	for(int m = 1; m <= n; m++)
	  {
	    Re_ei = exp_a_x_re[m];
	    Im_ei = exp_a_x_im[m];
	    Im_c = r_local_ci[index_start - m];
	    Re_c = r_local_ci[index_start+m];
	    r_local_ci[index_start - m] = Re_ei*Im_c - Im_ei*Re_c;
	    r_local_ci[index_start+m] = Re_ei*Re_c+Im_ei*Im_c;
	  }
      }

    //Step 7:       swap x and z
    swap_x_z(r_local_ci);

    //Step 8:       Final rotation around z - axis by(-alpha_z)
    for(int n = 1; n <= P; n++)
      {
	int index_start = n*n+n;
	real_t Re_ei, Im_ei, Im_c, Re_c;
	real_t oddoreven=-1;
	for(int m = 1; m <= n; m++)
	  {
	    Re_ei = exp_a_z_re[m];
	    Im_ei = exp_a_z_im[m];
	    Im_c = r_local_ci[index_start - m];
	    Re_c = r_local_ci[index_start+m];
	    r_local_ci[index_start - m] = Re_ei*Im_c - Im_ei*Re_c;
	    r_local_ci[index_start+m] = Re_ei*Re_c+Im_ei*Im_c;
	  }
      }

    //Now local expansion is expressed consistently with the original coordinates
    // Store the results in real - valued form

    for(int indice = 0; indice < NTERM; indice++)
      {
	#pragma omp atomic
	Ci->L[indice] += r_local_ci[indice];
      }
  }

  void M2L(Cell*Ci, Cell*Cj)
  {
    if(Ci->has_sink)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - Cj->X[d];

	real_t r_multipole_Cj[NTERM], r_local[NTERM];
	complex_t c_Cj[NTERM], c_local[NTERM];

	for(int index = 0; index < NTERM; index++)
	  {
	    c_local[index] = 0;
	    r_multipole_Cj[index] = Cj->M[index];
	  }

	for(int n = 1; n <= P; n++)
	  {
	    int index_start = n*n+n;
	    for(int m = 1; m <= n; m++)
	      r_multipole_Cj[index_start - m] *= -1;
	  }

	real_2_complex(r_multipole_Cj, c_Cj, P);

	std::vector < complex_t > Tnm(NTERM, 0);

	make_Tnm(&dX[0], &Tnm[0], P);

	for(int n = 0; n <= P; n++)
	  {
	    for(int m = -n; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= P - n; k++)
		  {
		    for(int l = -k; l <= k; l++)
		      {
			int indice2 = k*k+k+l;
			c_local[indice] += c_Cj[indice2]*Tnm[(n+k)*(n+k)+(n+k)+(m+l)];
		      }
		  }
	      }
	  }

	complex_2_real(c_local, r_local, P);

	for(int indice = 0; indice < NTERM; indice++)
	  {
#pragma omp atomic
	    Ci->L[indice] += r_local[indice];
	  }
      }
  }



  //below are the kernels for low order estimate of scalar forces
  //notice in P2M and M2M, all particles are sources and targets
  //such that the coefficients for the cells will be momentum-conserving

  void P2M_low(Cell*C)
  {
    complex_t c_multipole[9];
    for(Body*B = C->BODY; B != C->BODY+C->NBODY; B++)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = B->X[d] - C->X[d];
	complex_t Gnm[9];
	make_Gnm(&dX[0], &Gnm[0], 2);
	for(int indice = 0; indice < 9; indice++)
	  c_multipole[indice] += B->q*Gnm[indice];
      }
    real_t r_multipole[9];
    complex_2_real(&c_multipole[0], &r_multipole[0], 2);
    for(int indice = 0; indice < 9; indice++)
      C->M[indice] += r_multipole[indice];
  }

  void M2M_low(Cell*Ci)
  {
    complex_t c_multipole[9];
    for(int i = 0; i < 9; i++)
      {
	c_multipole[i] = 0;
      }

    for(Cell*Cj = Ci->CHILD; Cj != Ci->CHILD+Ci->NCHILD; Cj++)
      {
	for(int d = 0; d < 3; d++)
	  {
	    dX[d] = Cj->X[d] - Ci->X[d];
	  }

	real_t r_multipole_Cj[9];
	complex_t c_multipole_Cj[9];

	for(int i = 0; i < 9; i++)
	  {
	    r_multipole_Cj[i] = Cj->M[i];
	    c_multipole_Cj[i] = 0;
	  }

	real_2_complex(r_multipole_Cj, c_multipole_Cj, 2);

	complex_t Gnm[9];
	make_Gnm(&dX[0], &Gnm[0], 2);

	for(int n = 0; n <= 2; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= n; k++)
		  for(int l = std::max(-k, m - n+k); l <= std::min(k, m+n - k); l++)
		    c_multipole[indice] +=
		      c_multipole_Cj[(n - k)*(n - k)+(n - k)+(m - l)]*Gnm[k*k+k+l];
	      }
	  }
      }

    real_t r_multipole[9];
    complex_2_real(c_multipole, r_multipole, 2);

    for(int indice = 0; indice < 9; indice++)
      Ci->M[indice] += r_multipole[indice];
  }


  void M2L_low(Cell*Ci, Cell*Cj)
  {

    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];

    real_t r_multipole_Cj[9], r_local[9];
    complex_t c_Cj[9], c_local[9];

    for(int index = 0; index < 9; index++)
      {
	c_local[index] = 0;
	r_multipole_Cj[index] = Cj->M[index];
      }

    for(int n = 0; n <= 2; n++)
      {
	int index_start = n*n+n;
	for(int m = 1; m <= n; m++)
	  r_multipole_Cj[index_start - m] *= -1;
      }

    real_2_complex(r_multipole_Cj, c_Cj, 2);

    std::vector < complex_t > Tnm(9, 0);

    make_Tnm(&dX[0], &Tnm[0], 2);

    for(int n = 0; n <= 2; n++)
      {
	for(int m = 0; m <= n; m++)
	  {
	    int indice = n*n+n+m;
	    for(int k = 0; k <= 2 - n; k++)
	      {
		for(int l = -k; l <= k; l++)
		  {
		    int indice2 = k*k+k+l;
		    c_local[indice] += c_Cj[indice2]*Tnm[(n+k)*(n+k)+(n+k)+(m+l)];
		  }
	      }
	  }
      }

    complex_2_real(c_local, r_local, 2);

    for(int indice = 0; indice < 9; indice++)
      {
#pragma omp atomic
	Ci->L[indice] += r_local[indice];
      }
  }

  void P2P_low(Cell*Ci, Cell*Cj)
  {
    Body *Bi = Ci->BODY;
    Body *Bj = Cj->BODY;

    int ni = Ci->NBODY;
    int nj = Cj->NBODY;

#if SIMD_P2P
    int nii = (ni+NSIMD - 1) & (-NSIMD);

#ifndef DOUBLE_P2P
    float Xi[nii] __attribute__ ((aligned(64)));
    float Yi[nii] __attribute__ ((aligned(64)));
    float Zi[nii] __attribute__ ((aligned(64)));
#else
    double Xi[nii] __attribute__ ((aligned(64)));
    double Yi[nii] __attribute__ ((aligned(64)));
    double Zi[nii] __attribute__ ((aligned(64)));
#endif

    for(int k = 0; k < ni; k++)
      {
	Xi[k] = Bi[k].X[0];
	Yi[k] = Bi[k].X[1];
	Zi[k] = Bi[k].X[2];
      }

#ifndef DOUBLE_P2P
    Vec16f xi, yi, zi, r, mj, factor1, dx, dy, dz, r2;
#else
    Vec8d xi, yi, zi, r, mj, factor1, dx, dy, dz, r2;
#endif

    for(int i = 0; i < nii; i = i+NSIMD)
      {
	xi.load(Xi+i);
	yi.load(Yi+i);
	zi.load(Zi+i);
	for(int j = 0; j < nj; j++)
	  {
	    mj = Bj[j].q;
	    dx = xi - Bj[j].X[0];
	    dy = yi - Bj[j].X[1];
	    dz = zi - Bj[j].X[2];
	    r2 = dx*dx+dy*dy+dz*dz;
	    mj = select(r2 > 0, mj, 0);
	    r2 = select(r2 > 0, r2, HUGE);
	    factor1 = mj / r2;
	  }

	for(int k = 0; (k < NSIMD) && (i+k < ni); k++)
	  {
#pragma omp atomic
	    Bi[i+k].acc_old += (real_t) factor1[k];
	  }
      }
#else
    for(int i = 0; i < ni; i++)
      {
	real_t invR2 = 0;

	for(int j = 0; j < nj; j++)
	  {
	    for(int d = 0; d < 3; d++)
	      dX[d] = Bi[i].X[d] - Bj[j].X[d];

	    real_t R2 = norm(dX);

	    if(R2 > 0)
	      invR2 = Bj[j].q / R2;
	  }

#pragma omp atomic
	Bi[i].acc_old += invR2;
      }
#endif
  }

  void L2L_low(Cell*Ci)
  {
    real_t r_local_Ci[9];
    complex_t c_local_Ci[9];

    for(int indice = 0; indice < 9; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, 2);

    for(Cell*Cj = Ci->CHILD; Cj != Ci->CHILD+Ci->NCHILD; Cj++)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - Cj->X[d];

	complex_t c_local[9];

	for(int indice = 0; indice < 9; indice++)
	  c_local[indice] = 0.0;

	complex_t Gnm[9];

	make_Gnm_conj(&dX[0], &Gnm[0], 2);

	for(int n = 0; n <= 2; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= 2 - n; k++)
		  for(int l = -k; l <= k; l++)
		    c_local[indice] += c_local_Ci[(n+k)*(n+k)+(n+k)+(m+l)]*Gnm[k*k+k+l];
	      }
	  }

	real_t r_local[9];

	complex_2_real(c_local, r_local, 2);

	for(int indice = 0; indice < 9; indice++)
	  {
#pragma omp atomic
	    Cj->L[indice] += r_local[indice];
	  }

      }
  }


  void L2P_low(Cell*Ci)
  {
    real_t r_local_Ci[9];
    complex_t c_local_Ci[9];

    for(int indice = 0; indice < 9; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, 2);

    for(Body*B = Ci->BODY; B != Ci->BODY+Ci->NBODY; B++)
      {
	complex_t Phi[4];

	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - B->X[d];

	complex_t Gnm[9];
          
	make_Gnm_conj(&dX[0], &Gnm[0], 2);

	for(int n = 0; n <= 2; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n*n+n+m;
		for(int k = 0; k <= 2 - n; k++)
		  for(int l = -k; l <= k; l++)
		    Phi[indice] += c_local_Ci[(n+k)*(n+k)+(n+k)+(m+l)]*Gnm[k*k+k+l];
	      }
	  }

	B->acc_old += sqrt(std::real(Phi[3])*std::real(Phi[3]) +
                       std::imag(Phi[3])*std::imag(Phi[3]) +
                       std::real(Phi[2])*std::real(Phi[2]));
      }
  }
}
