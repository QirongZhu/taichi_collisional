#ifndef kernel_best_lean_h
#define kernel_best_lean_h

#include "exafmm.h"
#include <iostream>
#include <math.h>		/*isnan, sqrt, atan2 */
#include <stdio.h>		/*printf */

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

double dt_param = 0.025;

namespace exafmm
{
  real_t t_design[64 * 3] = {
    -.3345022001069780, .2559037493280540, -.9069848671303380,
    .2559037493281330, -.9069848671303220, -.3345022001069590,
    .3345022001069660, -.2559037493280760, -.9069848671303360,
    -.9069848671303240, -.3345022001070400, .2559037493280210,
    -.2559037493280450, -.9069848671303130, .3345022001070510,
    -.9069848671303130, .3345022001069890, -.2559037493281280,
    -.2559037493280770, .9069848671303280, -.3345022001069860,
    .9069848671303340, .3345022001069990, .2559037493280370,
    .2559037493280530, .9069848671303270, .3345022001070080,
    .9069848671303220, -.3345022001070060, -.2559037493280710,
    -.3345022001070360, -.2559037493280360, .9069848671303210,
    .3345022001069960, .2559037493281340, .9069848671303080,
    .1812453816504250, .1532808703053070, -.9714191095654030,
    .1532808703053200, -.9714191095653860, .1812453816505030,
    -.1812453816503960, -.1532808703053140, -.9714191095654070,
    -.9714191095654200, .1812453816504160, .1532808703052080,
    -.1532808703052110, -.9714191095654300, -.1812453816503590,
    -.9714191095653910, -.1812453816504850, -.1532808703053120,
    -.1532808703053150, .9714191095654040, .1812453816504130,
    .9714191095654130, -.1812453816504040, .1532808703052690,
    .1532808703052900, .9714191095654090, -.1812453816504040,
    .9714191095653960, .1812453816504510, -.1532808703053180,
    .1812453816504140, -.1532808703051840, .9714191095654240,
    -.1812453816504980, .1532808703053040, .9714191095653890,
    -.5671272467360530, .2239533904880130, .7925979844148650,
    .2239533904880520, .7925979844149510, -.5671272467359180,
    .5671272467359560, -.2239533904879880, .7925979844149420,
    .7925979844149380, -.5671272467359511, .2239533904880140,
    -.2239533904880400, .7925979844149330, .5671272467359481,
    .7925979844149070, .5671272467359810, -.2239533904880470,
    -.2239533904880380, -.7925979844149900, -.5671272467358700,
    -.7925979844149520, .5671272467359360, .2239533904880050,
    .2239533904880500, -.7925979844148791, .5671272467360200,
    -.7925979844148781, -.5671272467360270, -.2239533904880330,
    -.5671272467359670, -.2239533904880070, -.7925979844149290,
    .5671272467359780, .2239533904879970, -.7925979844149240,
    -0.6608756503404201, -.5307337939386410, -.8449580307443481,
    -.5307337939385640, -.8449580307444050, -0.6608756503393801,
    0.6608756503400100, .5307337939386150, -.8449580307443680,
    -.8449580307443810, -0.6608756503398400, -.5307337939385970,
    .5307337939386250, -.8449580307443620, 0.6608756503400000,
    -.8449580307444290, 0.6608756503394300, .5307337939385250,
    .5307337939386140, .8449580307443740, -0.6608756503393801,
    .8449580307443800, 0.6608756503393200, -.5307337939386040,
    -.5307337939386390, .8449580307443550, 0.6608756503397401,
    .8449580307443840, -0.6608756503398600, .5307337939385910,
    -0.6608756503402601, .5307337939385840, .8449580307443860,
    0.6608756503400399, -.5307337939385089, .8449580307444350,
    .6566234254592690, -.3208553693050780, -.6825668532284860,
    -.3208553693049820, -.6825668532284280, .6566234254593770,
    -.6566234254592740, .3208553693050560, -.6825668532284920,
    -.6825668532284440, .6566234254593130, -.3208553693050770,
    .3208553693049670, -.6825668532284750, -.6566234254593351,
    -.6825668532284180, -.6566234254593410, .3208553693050760,
    .3208553693050100, .6825668532284710, .6566234254593180,
    .6825668532284540, -.6566234254593000, -.3208553693050810,
    -.3208553693049920, .6825668532284720, -.6566234254593259,
    .6825668532284560, .6566234254593050, .3208553693050690,
    .6566234254593371, .3208553693050400, .6825668532284380,
    -.6566234254593390, -.3208553693050600, .6825668532284270,
    -.5773502691896160, -.5773502691896200, -.5773502691896411,
    .5773502691896250, .5773502691896060, -.5773502691896460,
    .5773502691896440, -.5773502691895930, .5773502691896401,
    -.5773502691896371, .5773502691896250, .5773502691896160
  };

  real_t oddOrEventable[EXPANSION + 1];

  real_t factorial_table[30], inverse_factorial_table[30];

  real_t factorial_coef[EXPANSION + 1][EXPANSION + 1];
  real_t factorial_coef_inv[EXPANSION + 1][EXPANSION + 1];

  real_t factorial_coef_oned[(EXPANSION + 1) * (EXPANSION + 1)];
  real_t factorial_coef_inv_oned[(EXPANSION + 1) * (EXPANSION + 1)];

  real_t combinator_coef[EXPANSION + 1][EXPANSION + 1];

  const complex_t R(1., 0.);	//!<Real unit
  const complex_t I(0., 1.);	//!<Imaginary unit

  real_t norm(real_t * X)
  {
    return X[0] * X[0] + X[1] * X[1] + X[2] * X[2];
  }


  int oddOrEven(int n)
  {
    return (((n) & 1) == 1) ? -1 : 1;
  }

  void initKernel()
  {
    NTERM = (P + 1) * (P + 1);
    factorial_table[0] = 1.0;
    inverse_factorial_table[0] = 1.0;
    for(int i = 1; i < 30; i++)
      {
	factorial_table[i] = (real_t) i *factorial_table[i - 1];
	inverse_factorial_table[i] = 1 / factorial_table[i];
      }
    for(int i = 0; i <= P; i++)
      {
	for(int j = 0; j <= P; j++)
	  {
	    factorial_coef[i][j] = factorial_table[i + j] * factorial_table[i - j];
	    factorial_coef_inv[i][j] = 1.0 / factorial_coef[i][j];
	    if(j <= i)
	      combinator_coef[i][j] =
		factorial_table[i] * inverse_factorial_table[j] * inverse_factorial_table[i - j];
	  }

	oddOrEventable[i] = oddOrEven(i);

      }

    for(int n = 0; n <= P; n++)
      {
	int index_start = n * n + n;
	factorial_coef_oned[index_start] = factorial_coef[n][0];
	factorial_coef_inv_oned[index_start] = factorial_coef_inv[n][0];
	for(int m = 1; m <= n; m++)
	  {
	    factorial_coef_oned[index_start + m] = factorial_coef[n][m];
	    factorial_coef_oned[index_start - m] = factorial_coef[n][m];
	    factorial_coef_inv_oned[index_start + m] = factorial_coef_inv[n][m];
	    factorial_coef_inv_oned[index_start - m] = factorial_coef_inv[n][m];
	  }
      }
  }

  void real_2_complex(real_t * real_arr, complex_t * complex_arr, int order)
  {
    for(int n = 0; n <= order; n++)
      {
	int index_start = n * n + n;
	complex_arr[index_start] = complex_t(real_arr[index_start], 0);
	real_t oddevenfac = -1;
	for(int m = 1; m <= n; m++)
	  {
	    complex_arr[index_start - m] = complex_t(real_arr[index_start + m] * oddevenfac,
						     -real_arr[index_start - m] * oddevenfac);
	    complex_arr[index_start + m] = complex_t(real_arr[index_start + m], real_arr[index_start - m]);
	    oddevenfac *= -1;
	  }
      }
  }

  void complex_2_real(complex_t * complex_arr, real_t * real_arr, int order)
  {
    for(int n = 0; n <= order; n++)
      {
	int index_start = n * n + n;
	real_arr[index_start] = std::real(complex_arr[index_start]);
	for(int m = 1; m <= n; m++)
	  {
	    real_arr[index_start - m] = std::imag(complex_arr[index_start + m]);
	    real_arr[index_start + m] = std::real(complex_arr[index_start + m]);
	  }
      }
  }

  void make_Tnm(real_t * dX, complex_t * Tnm, int Order)
  {
    real_t r2 = norm(dX);
    real_t invr2 = 1 / r2;
    Tnm[0] = complex_t(sqrt(invr2), 0);

    for(int n = 1; n <= Order; n++)
      {
	int n2 = n * n;
	int ind = n2 + n + n;

	Tnm[ind] = Tnm[n2 - 1] * complex_t(dX[0] * (2 * n - 1) * invr2, dX[1] * (2 * n - 1) * invr2);

	if(n > 0)
	  {
	    int m, m2, indexstart = n2 - n, indexnew = n2 + n;

	    if(n > 1)
	      {
		for(m = 0; m < n - 1; m++)
		  {
		    m2 = m * m;
		    Tnm[indexnew + m] =
		      invr2 * (2 * n - 1.) * dX[2] * Tnm[indexstart + m] -
		      invr2 * ((n - 1) * (n - 1) * (n > 1) - m * m) * Tnm[n2 - 3 * n + 2 + m];
		  }
	      }

	    m = n - 1;
	    Tnm[indexnew + m] = invr2 * (2 * n - 1.) * dX[2] * Tnm[indexstart + m];
	  }

	real_t oddoreven = -1;

	for(int m = 1; m <= n; m++)
	  {
	    Tnm[n2 + n - m] = std::conj(Tnm[n2 + n + m] * oddoreven);
	    oddoreven *= -1;
	  }
      }
  }

  void make_Gnm(real_t * dX, complex_t * Gnm, int Order)
  {
    real_t r2 = norm(dX);

    Gnm[0] = complex_t(1, 0);

    for(int n = 1; n <= Order; n++)
      {
	int n2 = n * n;
	int ind = n2 + n + n;

	Gnm[ind] = Gnm[n2 - 1] * complex_t(dX[0] / (2 * n), dX[1] / (2 * n));

	if(n > 0)
	  {
	    int m, m2, indexstart = n2 - n, indexnew = n2 + n;

	    if(n > 1)
	      {
		for(m = 0; m < n - 1; m++)
		  {
		    m2 = m * m;
		    Gnm[indexnew + m] =
		      (2 * n - 1.) / (n2 - m2) * dX[2] * Gnm[indexstart + m] - r2 / (n2 - m2) * Gnm[n2 -
												    3 * n +
												    2 + m];
		  }
	      }

	    m = n - 1;
	    Gnm[indexnew + m] = dX[2] * Gnm[indexstart + m];
	  }

	real_t oddoreven = -1;
	for(int m = 1; m <= n; m++)
	  {
	    Gnm[n2 + n - m] = std::conj(Gnm[n2 + n + m] * oddoreven);
	    oddoreven *= -1;
	  }
      }
  }

  void P2P_simple(Cell * Ci, Cell * Cj)
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
	    int nii = ceil((float) ni / (float) NSIMD) * NSIMD;

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
		Xi[k] = -Bi[k].X[0];
		Yi[k] = -Bi[k].X[1];
		Zi[k] = -Bi[k].X[2];
	      }

#ifndef DOUBLE_P2P
	    Vec16f xi, yi, zi, r, mj;
	    Vec16f invR, factor1, fac1, dx, dy, dz, r2;
	    Vec16f ax, ay, az, pot;
#else
	    Vec8d xi, yi, zi, r, mj;
            Vec8d invR, factor1, fac1, dx, dy, dz, r2;
            Vec8d ax, ay, az, pot;
#endif
	    for(int i = 0; i < nii; i = i + NSIMD)
	      {
		xi.load(Xi + i);
		yi.load(Yi + i);
		zi.load(Zi + i);

		ax = 0, ay = 0, az = 0, pot = 0;

		for(int j = 0; j < nj; j++)
		  {
		    if(!Bj[j].issource)
		      continue;

		    mj = Bj[j].q;
		    dx = Bj[j].X[0] + xi;
		    dy = Bj[j].X[1] + yi;
		    dz = Bj[j].X[2] + zi;

		    r2 = dx * dx + dy * dy + dz * dz;

		    mj = select(r2 > 0, mj, 0);
		    r2 = select(r2 > 0, r2, 1e38);

		    r = sqrt(r2);
		    invR = 1 / r;
		    mj *= invR;

		    pot += mj;
		    factor1 = mj * (invR * invR);
		    ax += dx * factor1;
		    ay += dy * factor1;
		    az += dz * factor1;
		  }

		for(int k = 0; (k < NSIMD) && (i + k < ni); k++)
		  {
		    if(Bi[i + k].issink)
		      {
#pragma omp atomic
			Bi[i + k].p += (real_t) pot[k];
#pragma omp atomic
			Bi[i + k].F[0] += (real_t) ax[k];
#pragma omp atomic
			Bi[i + k].F[1] += (real_t) ay[k];
#pragma omp atomic
			Bi[i + k].F[2] += (real_t) az[k];
		      }
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

		real_t ax = 0;
		real_t ay = 0;
		real_t az = 0;
		real_t pot = 0;

		for(int j = 0; j < nj; j++)
		  {
		    if(!Bj[j].issource)
		      continue;

		    for(int d = 0; d < 3; d++)
		      {
			dX[d] = Bj[j].X[d] - Bi[i].X[d];
		      }

		    real_t R2 = norm(dX);

		    if(R2 > 0)
		      {
			real_t R = sqrt(R2);

			real_t invR2 = 1 / R2;
			real_t invR = Bj[j].q * sqrt(invR2) * Bj[j].issource;

			pot += invR;

			for(int d = 0; d < 3; d++)
			  dX[d] *= invR2 * invR;
			ax += dX[0];
			ay += dX[1];
			az += dX[2];
		      }
		  }

#pragma omp atomic
		Bi[i].p += (real_t) pot;
#pragma omp atomic
		Bi[i].F[0] += (real_t) ax;
#pragma omp atomic
		Bi[i].F[1] += (real_t) ay;
#pragma omp atomic
		Bi[i].F[2] += (real_t) az;
	      }
	  }
#else
	// do not use vectorized version if n1*n2 too small
	for(int i = 0; i < ni; i++)
	  {
	    real_t ax = 0;
	    real_t ay = 0;
	    real_t az = 0;
	    real_t pot = 0;

	    for(int j = 0; j < nj; j++)
	      {
		for(int d = 0; d < 3; d++)
		  dX[d] = Bi[i].X[d] - Bj[j].X[d];

		real_t R2 = norm(dX);

		if(R2 > 0)
		  {
		    real_t R = sqrt(R2);

		    real_t invR2 = 1.0 / R2;
		    real_t invR = Bj[j].q * sqrt(invR2) * Bj[j].issource;

		    pot += invR;

		    for(int d = 0; d < 3; d++)
		      dX[d] *= invR2 * invR;

		    ax += dX[0];
		    ay += dX[1];
		    az += dX[2];
		  }
	      }

	    if(Bi[i].issink)
	      {
#pragma omp atomic
		Bi[i].p += pot;
#pragma omp atomic
		Bi[i].F[0] -= ax;
#pragma omp atomic
		Bi[i].F[1] -= ay;
#pragma omp atomic
		Bi[i].F[2] -= az;
	      }
	  }
#endif


#pragma omp atomic
	Ci->NP2P += 1;

      }
  }


  void P2P(Cell * Ci, Cell * Cj)
  {
    if(Ci->has_sink)
      {
	Body *Bi = Ci->BODY;
	Body *Bj = Cj->BODY;

	int ni = Ci->NBODY;
	int nj = Cj->NBODY;

	real_t delta_t;

#if SIMD_P2P
	int nii = ceil((float) ni / (float) NSIMD) * NSIMD;

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
	Vec16f ax, ay, az, pot, timestep, tau, dtau;
#else
	Vec8d mi, xi, yi, zi, vxi, vyi, vzi, r2, v2, vdotdr2, r, mj;
        Vec8d invR, factor1, fac1, dx, dy, dz, dvx, dvy, dvz;
        Vec8d ax, ay, az, pot, timestep, tau, dtau;
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

	    ax = 0, ay = 0, az = 0, pot = 0, timestep = 1e38;

	    for(int j = 0; j < nj; j++)
	      {
		mj = Bj[j].q;
		dx = xi + Bj[j].X[0];
		dy = yi + Bj[j].X[1];
		dz = zi + Bj[j].X[2];
		dvx = vxi + Bj[j].V[0];
		dvy = vyi + Bj[j].V[1];
		dvz = vzi + Bj[j].V[2];

		r2 = dx * dx + dy * dy + dz * dz;
		v2 = dvx * dvx + dvy * dvy + dvz * dvz + 1e-20;

		mj = select(r2 > 0, mj, 0);
		mi += mj;
		r2 = select(r2 > 0, r2, 1e38);

		vdotdr2 = (dx * dvx + dy * dvy + dz * dvz) / r2;

		r = sqrt(r2);
		invR = 1 / r;

		tau = dt_param / M_SQRT2 * r * sqrt(r / mi);
		dtau = 3 * tau * vdotdr2 / 2;
		dtau = select(dtau < 1, dtau, 1);
		tau /= (1 - dtau / 2);
		timestep = min(tau, timestep);

		tau = dt_param * r / sqrt(v2);
		dtau = tau * vdotdr2 * (1 + mi / (v2 * r));
		dtau = select(dtau < 1, dtau, 1);
		tau /= (1 - dtau / 2);
		timestep = min(tau, timestep);

		mj *= invR;

		if(Bj[j].issource)
		  {
		    pot += mj;
		    factor1 = mj * (invR * invR);
		    ax += dx * factor1;
		    ay += dy * factor1;
		    az += dz * factor1;
		  }
	      }

	    for(int k = 0; (k < NSIMD) && (i + k < ni); k++)
	      {
		if(Bi[i + k].issink)
		  {
#pragma omp atomic
		    Bi[i + k].p    += (real_t) pot[k];
#pragma omp atomic
		    Bi[i + k].F[0] += (real_t) ax[k];
#pragma omp atomic
		    Bi[i + k].F[1] += (real_t) ay[k];
#pragma omp atomic
		    Bi[i + k].F[2] += (real_t) az[k];
		   
		    if(Bi[i + k].timestep > (real_t)timestep[k])
		      Bi[i + k].timestep = (real_t)timestep[k];	      
	     
		  }
	      }
	  }
#else
	for(int i = 0; i < ni; i++)
	  {
	    real_t ax = 0;
	    real_t ay = 0;
	    real_t az = 0;
	    real_t pot = 0;

	    real_t timestep = 1e38;

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
		    real_t R = sqrt(R2);
		    real_t vdotdr2;
		    real_t v2 = norm(dV);

		    vdotdr2 = (dX[0] * dV[0] + dX[1] * dV[1] + dX[2] * dV[2]) / R2;

		    real_t tau = dt_param / M_SQRT2 * sqrt(R * R2 / (Bi[i].q + Bj[j].q));
		    real_t dtau = 3 * tau * vdotdr2 / 2;
		    if(dtau > 1)
		      dtau = 1;
		    tau /= (1 - dtau / 2);
		    if(tau < timestep)
		      timestep = tau;

		    if(v2 > 0)
		      {
			tau = dt_param * R / sqrt(v2);
			dtau = tau * vdotdr2 * (1 + (Bi[i].q + Bj[j].q) / (v2 * R));
			if(dtau > 1)
			  dtau = 1;
			tau /= (1 - dtau / 2);
			if(tau < timestep)
			  timestep = tau;
		      }

		    real_t invR2 = 1.0 / R2;
		    real_t invR = Bj[j].q * sqrt(invR2) * Bj[j].issource;

		    pot += invR;

		    for(int d = 0; d < 3; d++)
		      dX[d] *= invR2 * invR;

		    ax += dX[0];
		    ay += dX[1];
		    az += dX[2];
		  }
	      }

	    if(Bi[i].issink)
	      {
#pragma omp atomic
		Bi[i].p += pot;
#pragma omp atomic
		Bi[i].F[0] += ax;
#pragma omp atomic
		Bi[i].F[1] += ay;
#pragma omp atomic
		Bi[i].F[2] += az;

		if(Bi[i].timestep > timestep)
		  Bi[i].timestep = timestep;

	      }
	  }
#endif

#pragma omp atomic
	Ci->NP2P += 1;
      }
  }

  void P2M(Cell * C)
  {
    real_t min_acc = 1e30;

#ifndef MINIBALL
    if(C->NBODY > 0)
      {
	real_t comx = 0, comy = 0, comz = 0, totalq = 0;

	for(Body * B = C->BODY; B != C->BODY + C->NBODY; B++)
	  {
	    comx += B->X[0] * B->q;
	    comy += B->X[1] * B->q;
	    comz += B->X[2] * B->q;
	    totalq += B->q;

	    min_acc = fmin(B->acc_old, min_acc);
	    if(B->issink)
	      C->has_sink = true;
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
	C->R = 0;
      }
#else
    if(C->NBODY > 0)
      {
	std::list < std::vector < real_t > >lp;

	for(Body * B = C->BODY; B != C->BODY + C->NBODY; B++)
	  {
	    std::vector < real_t > p(3);
	    p[0] = B->X[0];
	    p[1] = B->X[1];
	    p[2] = B->X[2];
	    lp.push_back(p);

	    min_acc = fmin(B->acc_old, min_acc);

	    if(B->issink)
	      C->has_sink = true;
	  }
	// define the types of iterators through the points and their coordinates
	typedef std::list < std::vector < real_t > >::const_iterator PointIterator;
	typedef std::vector < real_t >::const_iterator CoordIterator;
	typedef Miniball:: Miniball < Miniball::CoordAccessor < PointIterator, CoordIterator > > MB;
	MB mb(3, lp.begin(), lp.end());
	real_t rminball = sqrt(mb.squared_radius());
	const real_t *center = mb.center();

	if(rminball < C->R)
	  {			//use rminball as R
	    C->X[0] = center[0];
	    C->X[1] = center[1];
	    C->X[2] = center[2];
	    C->R = rminball + 1e-6;
	  }
      }
    else
      {
	C->R = 0;
      }  
#endif

    complex_t c_multipole[NTERM];

    for(Body * B = C->BODY; B != C->BODY + C->NBODY; B++)
      {
	for(int d = 0; d < 3; d++)
	  {
	    dX[d] = B->X[d] - C->X[d];
	  }

	complex_t Gnm[NTERM];

	make_Gnm(&dX[0], &Gnm[0], P);

	for(int indice = 0; indice < NTERM; indice++)
	  c_multipole[indice] += (B->q * B->issource) * Gnm[indice];
      }

    real_t r_multipole[NTERM];
    complex_2_real(&c_multipole[0], &r_multipole[0], P);

    for(int indice = 0; indice < NTERM; indice++)
      C->M[indice] += r_multipole[indice];

    C->min_acc = min_acc;
  }

  void M2M(Cell * Ci)
  {
    real_t min_acc = 1e30;

#ifndef MINIBALL
    if(Ci->NCHILD > 0)
      {
	real_t comx = 0, comy = 0, comz = 0, totalq = 0;

	for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
	  {
	    comx += Cj->X[0] * Cj->cell_mass;
	    comy += Cj->X[1] * Cj->cell_mass;
	    comz += Cj->X[2] * Cj->cell_mass;
	    totalq += Cj->cell_mass;
	    min_acc = fmin(Cj->min_acc, min_acc);

	    if(Cj->has_sink)
	      Ci->has_sink = true;
	  }

	comx /= totalq;
	comy /= totalq;
	comz /= totalq;

	Ci->X[0] = comx;
	Ci->X[1] = comy;
	Ci->X[2] = comz;
	Ci->cell_mass = totalq;

	real_t max_r = 1e-10;

	for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
	  {
	    for(int d = 0; d < 3; d++)
	      dX[d] = Cj->X[d] - Ci->X[d];
	    max_r = fmax(max_r, sqrt(norm(dX)) + Cj->R);
	  }

	Ci->R = max_r;
      }
    else
      {
	Ci->R = 0;
      }
#else
    if (Ci->NCHILD > 0) {

      std::list<std::vector<real_t> > lp;

      for (Cell *ci=Ci->CHILD; ci!=Ci->CHILD+Ci->NCHILD; ci++) {

	min_acc = fmin(ci->min_acc, min_acc);

	if(ci->has_sink) Ci->has_sink=true;

	if(ci->NCHILD == 0) { //leaves, just use bodies to get R_max
	  for(Body*B=ci->BODY; B!=ci->BODY+ci->NBODY; B++) {
	    std::vector<real_t> p(3);
	    p[0] = B->X[0];
	    p[1] = B->X[1];
	    p[2] = B->X[2];
	    lp.push_back(p);
	  }
	}
	else{ 
	  // use granddaughers information to get R_max with a 
	  // precomputed t-design of a sphere with 64 points.
	  for (Cell *cii=ci->CHILD; cii!=ci->CHILD+ci->NCHILD; cii++) {
	    for(int k=0; k<64; k++) {
	      std::vector<real_t> p(3);
	      p[0] = cii->X[0] + cii->R * t_design[3*k+0];
	      p[1] = cii->X[1] + cii->R * t_design[3*k+1];
	      p[2] = cii->X[2] + cii->R * t_design[3*k+2];
	      lp.push_back(p);
	    }
	  }
	}
      }

      typedef std::list<std::vector<real_t> >::const_iterator PointIterator;
      typedef std::vector<real_t>::const_iterator CoordIterator;
      typedef Miniball::
	Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > MB;
      MB mb (3, lp.begin(), lp.end());

      real_t rminball = sqrt(mb.squared_radius());

      if(Ci->R  > rminball) {
	const real_t* center = mb.center();
	Ci->X[0]    = center[0];
	Ci->X[1]    = center[1];
	Ci->X[2]    = center[2];
	Ci->R       = rminball + 1e-6;
      }
    
      Ci->min_acc = min_acc;
    }
    else
      {
	Ci->R  = 0;
      }
#endif

    complex_t c_multipole[NTERM];
    for(int i = 0; i < NTERM; i++)
      {
	c_multipole[i] = 0;
      }

    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
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
		int indice = n * n + n + m;
		for(int k = 0; k <= n; k++)
		  for(int l = std::max(-k, m - n + k); l <= std::min(k, m + n - k); l++)
		    c_multipole[indice] +=
		      c_multipole_Cj[(n - k) * (n - k) + (n - k) + (m - l)] * Gnm[k * k + k + l];
	      }
	  }
      }

    real_t r_multipole[NTERM];
    complex_2_real(c_multipole, r_multipole, P);

    for(int indice = 0; indice < NTERM; indice++)
      Ci->M[indice] += r_multipole[indice];

    Ci->min_acc = min_acc;
  }

  void L2L(Cell * Ci)
  {
    real_t r_local_Ci[NTERM];
    complex_t c_local_Ci[NTERM];

    for(int indice = 0; indice < NTERM; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, P);

    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - Cj->X[d];

	complex_t c_local[NTERM];

	for(int indice = 0; indice < NTERM; indice++)
	  c_local[indice] = 0.0;

	complex_t Gnm[NTERM];

	make_Gnm(&dX[0], &Gnm[0], P);

	for(int i = 0; i < NTERM; i++)
	  Gnm[i] = std::conj(Gnm[i]);

	for(int n = 0; n <= P; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n * n + n + m;
		for(int k = 0; k <= P - n; k++)
		  for(int l = -k; l <= k; l++)
		    c_local[indice] += c_local_Ci[(n + k) * (n + k) + (n + k) + (m + l)] * Gnm[k * k + k + l];
	      }
	  }

	real_t r_local[NTERM];

	complex_2_real(c_local, r_local, P);

	for(int indice = 0; indice < NTERM; indice++)
	  Cj->L[indice] += r_local[indice];

      }
  }

  void L2P(Cell * Ci)
  {
    real_t r_local_Ci[NTERM];
    complex_t c_local_Ci[NTERM];

    for(int indice = 0; indice < NTERM; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, P);

    for(Body * B = Ci->BODY; B != Ci->BODY + Ci->NBODY; B++)
      {
	complex_t Phi[4];

	if(!B->issink)
	  continue;

	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - B->X[d];

	complex_t Gnm[NTERM];
	make_Gnm(&dX[0], &Gnm[0], P);

	for(int i = 0; i < NTERM; i++)
	  Gnm[i] = std::conj(Gnm[i]);

	for(int n = 0; n <= 1; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n * n + n + m;

		if(indice == 1)
		  continue;

		for(int k = 0; k <= P - n; k++)
		  for(int l = -k; l <= k; l++)
		    Phi[indice] += c_local_Ci[(n + k) * (n + k) + (n + k) + (m + l)] * Gnm[k * k + k + l];
	      }
	  }

	B->p += std::real(Phi[0]);
	B->F[0] -= std::real(Phi[3]);
	B->F[1] -= std::imag(Phi[3]);
	B->F[2] -= std::real(Phi[2]);
      }
  }

  void M2L_rotate(Cell * Ci, Cell * Cj)
  {
    if(Ci->has_sink)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - Cj->X[d];

	real_t r2_xy = dX[0] * dX[0] + dX[1] * dX[1];
	real_t r_xy = sqrt(r2_xy);
	real_t r2 = r2_xy + dX[2] * dX[2];
	real_t r = sqrt(r2);
	real_t r_inv = 1.0 / r;

#if 1
	real_t exp_a_z_re[P + 1], exp_a_z_im[P + 1];
	real_t exp_a_x_re[P + 1], exp_a_x_im[P + 1];

	exp_a_z_re[0] = 1;
	exp_a_z_im[0] = 0;
	exp_a_x_re[0] = 1;
	exp_a_x_im[0] = 0;

	exp_a_z_re[1] = 1;
	exp_a_z_im[1] = 0;
	exp_a_x_re[1] = dX[2] * r_inv;
	exp_a_x_im[1] = -r_xy * r_inv;

	if(r_xy > 0)
	  {
	    r_xy = 1.0 / r_xy;
	    exp_a_z_re[1] = dX[1] * r_xy;
	    exp_a_z_im[1] = dX[0] * r_xy;
	  }

	for(int m = 1; m < P; m++)
	  {
	    exp_a_z_re[m + 1] = exp_a_z_re[m] * exp_a_z_re[1] - exp_a_z_im[m] * exp_a_z_im[1];
	    exp_a_z_im[m + 1] = exp_a_z_re[m] * exp_a_z_im[1] + exp_a_z_im[m] * exp_a_z_re[1];
	    exp_a_x_re[m + 1] = exp_a_x_re[m] * exp_a_x_re[1] - exp_a_x_im[m] * exp_a_x_im[1];
	    exp_a_x_im[m + 1] = exp_a_x_re[m] * exp_a_x_im[1] + exp_a_x_im[m] * exp_a_x_re[1];
	  }
#endif

	real_t r_multipole_Cj[NTERM], r_local[NTERM];

	//Step 0:       copy multipoles, and make multipoles into homogenius polynomials
	for(int indice = 0; indice < NTERM; indice++)
	  {
	    r_local[indice] = 0;
	    r_multipole_Cj[indice] = Cj->M[indice] * factorial_coef_oned[indice];
	  }

#if 1
	//Step 1:       first rotate multipoles around the original z - axis with alpha_z
	for(int n = 1; n <= P; n++)
	  {
	    int index_start = n * n + n;
	    real_t Re_ei, Im_ei, Im_r, Re_r;
	    for(int m = 1; m <= n; m++)
	      {
		Re_ei = exp_a_z_re[m];
		Im_ei = exp_a_z_im[m];
		Im_r = r_multipole_Cj[index_start - m];
		Re_r = r_multipole_Cj[index_start + m];
		r_multipole_Cj[index_start - m] = Re_ei * Im_r + Im_ei * Re_r;
		r_multipole_Cj[index_start + m] = Re_ei * Re_r - Im_ei * Im_r;
	      }
	  }
#endif

	//Step 2:       swap x and z
	swap_x_z(r_multipole_Cj);

#if 1
	//step 3:       rotate around(new) z - axis by alpha_x
	for(int n = 1; n <= P; n++)
	  {
	    int index_start = n * n + n;
	    real_t Re_ei, Im_ei, Im_r, Re_r;
	    for(int m = 1; m <= n; m++)
	      {
		Re_ei = exp_a_x_re[m];
		Im_ei = exp_a_x_im[m];
		Im_r = r_multipole_Cj[index_start - m];
		Re_r = r_multipole_Cj[index_start + m];
		r_multipole_Cj[index_start - m] = Re_ei * Im_r + Im_ei * Re_r;
		r_multipole_Cj[index_start + m] = Re_ei * Re_r - Im_ei * Im_r;
	      }
	  }
#endif

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
	    int index_start = n * (n + 1);
	    powers_of_r = powers_of_r_n;

	    for(int k = 0; k <= P - n; k++)
	      {
		real_t fac = factorial_table[n + k] * powers_of_r;
		r_local[index_start] += r_multipole_Cj[k * (k + 1)] * fac;
		powers_of_r *= r_inv;
	      }

	    powers_of_r_m = r_inv;

	    real_t oddevenfac = -1;

	    for(int m = 1; m <= std::min(n, P - n); m++)
	      {
		powers_of_r = powers_of_r_n * powers_of_r_m;
		for(int k = m; k <= P - n; k++)
		  {
		    real_t fac = factorial_table[n + k] * powers_of_r * oddevenfac;
		    int index_start_k = k * (k + 1);
		    r_local[index_start - m] += r_multipole_Cj[index_start_k - m] * fac;
		    r_local[index_start + m] += r_multipole_Cj[index_start_k + m] * fac;
		    powers_of_r *= r_inv;
		  }
		oddevenfac *= -1;
		powers_of_r_m *= r_inv;
	      }
	    powers_of_r_n *= r_inv;
	  }

	//Step 5:       swap x and z
	swap_x_z(r_local);

#if 1
	//Step 6:       Rotate local expansion around z - axis by(-alpha_x)
	for(int n = 1; n <= P; n++)
	  {
	    int index_start = n * n + n;
	    real_t Re_ei, Im_ei, Im_c, Re_c;
	    for(int m = 1; m <= n; m++)
	      {
		Re_ei = exp_a_x_re[m];
		Im_ei = exp_a_x_im[m];
		Im_c = r_local[index_start - m];
		Re_c = r_local[index_start + m];
		r_local[index_start - m] = Re_ei * Im_c - Im_ei * Re_c;
		r_local[index_start + m] = Re_ei * Re_c + Im_ei * Im_c;
	      }
	  }
#endif

	//Step 7:       swap x and z
	swap_x_z(r_local);

#if 1
	//Step 8:       Final rotation around z - axis by(-alpha_z)
	for(int n = 1; n <= P; n++)
	  {
	    int index_start = n * n + n;
	    real_t Re_ei, Im_ei, Im_c, Re_c;
	    for(int m = 1; m <= n; m++)
	      {
		Re_ei = exp_a_z_re[m];
		Im_ei = exp_a_z_im[m];
		Im_c = r_local[index_start - m];
		Re_c = r_local[index_start + m];
		r_local[index_start - m] = Re_ei * Im_c - Im_ei * Re_c;
		r_local[index_start + m] = Re_ei * Re_c + Im_ei * Im_c;
	      }
	  }
#endif
	//Now local expansion is expressed consistently with the original coordinates

	// Store the results in real - valued form              
	for(int indice = 0; indice < NTERM; indice++)
	  {
#pragma omp atomic
	    Ci->L[indice] += r_local[indice];
	  }
#pragma omp atomic
	Ci->NM2L += 1;
      }
  }

  void M2L(Cell * Ci, Cell * Cj)
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
	    int index_start = n * n + n;
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
		int indice = n * n + n + m;
		for(int k = 0; k <= P - n; k++)
		  {
		    for(int l = -k; l <= k; l++)
		      {
			int indice2 = k * k + k + l;
			c_local[indice] += c_Cj[indice2] * Tnm[(n + k) * (n + k) + (n + k) + (m + l)];
		      }
		  }
	      }
	  }

	complex_2_real(c_local, r_local, P);

	for(int indice = 0; indice < NTERM; indice++) {
#pragma omp atomic
	  Ci->L[indice] += r_local[indice];
	}
      }
  }



  //below are the kernels for low order estimate of scalar forces
  //notice in P2M and M2M, all particles are sources and targets
  //such that the coefficients for the cells will be momentum-conserving

  void P2M_low(Cell * C)
  {
    complex_t c_multipole[9];
    for(Body * B = C->BODY; B != C->BODY + C->NBODY; B++)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = B->X[d] - C->X[d];
	complex_t Gnm[9];
	make_Gnm(&dX[0], &Gnm[0], 2);
	for(int indice = 0; indice < 9; indice++)
	  c_multipole[indice] += B->q * Gnm[indice];
      }
    real_t r_multipole[9];
    complex_2_real(&c_multipole[0], &r_multipole[0], 2);
    for(int indice = 0; indice < 9; indice++)
      C->M[indice] += r_multipole[indice];
  }

  void M2M_low(Cell * Ci)
  {
    complex_t c_multipole[9];
    for(int i = 0; i < 9; i++)
      {
	c_multipole[i] = 0;
      }

    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
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
		int indice = n * n + n + m;
		for(int k = 0; k <= n; k++)
		  for(int l = std::max(-k, m - n + k); l <= std::min(k, m + n - k); l++)
		    c_multipole[indice] +=
		      c_multipole_Cj[(n - k) * (n - k) + (n - k) + (m - l)] * Gnm[k * k + k + l];
	      }
	  }
      }

    real_t r_multipole[9];
    complex_2_real(c_multipole, r_multipole, 2);

    for(int indice = 0; indice < 9; indice++)
      Ci->M[indice] += r_multipole[indice];
  }


  void M2L_low(Cell * Ci, Cell * Cj)
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
	int index_start = n * n + n;
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
	    int indice = n * n + n + m;
	    for(int k = 0; k <= 2 - n; k++)
	      {
		for(int l = -k; l <= k; l++)
		  {
		    int indice2 = k * k + k + l;
		    c_local[indice] += c_Cj[indice2] * Tnm[(n + k) * (n + k) + (n + k) + (m + l)];
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

  void P2P_low(Cell * Ci, Cell * Cj)
  {
    Body *Bi = Ci->BODY;
    Body *Bj = Cj->BODY;

    int ni = Ci->NBODY;
    int nj = Cj->NBODY;

#if SIMD_P2P
    int nii = ceil((float) ni / (float) NSIMD) * NSIMD;

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

    for(int i = 0; i < nii; i = i + NSIMD)
      {
	xi.load(Xi + i);
	yi.load(Yi + i);
	zi.load(Zi + i);
	for(int j = 0; j < nj; j++)
	  {
	    mj = Bj[j].q;
	    dx = xi - Bj[j].X[0];
	    dy = yi - Bj[j].X[1];
	    dz = zi - Bj[j].X[2];
	    r2 = dx * dx + dy * dy + dz * dz;
	    mj = select(r2 > 0, mj, 0);
	    r2 = select(r2 > 0, r2, HUGE);
	    factor1 = mj / r2;
	  }

	for(int k = 0; (k < NSIMD) && (i + k < ni); k++)
	  {
#pragma omp atomic
	    Bi[i + k].acc_old += (real_t) factor1[k];
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

  void L2L_low(Cell * Ci)
  {
    real_t r_local_Ci[9];
    complex_t c_local_Ci[9];

    for(int indice = 0; indice < 9; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, 2);

    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
      {
	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - Cj->X[d];

	complex_t c_local[9];

	for(int indice = 0; indice < 9; indice++)
	  c_local[indice] = 0.0;

	complex_t Gnm[9];

	make_Gnm(&dX[0], &Gnm[0], 2);

	for(int i = 0; i < 9; i++)
	  Gnm[i] = std::conj(Gnm[i]);

	for(int n = 0; n <= 2; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n * n + n + m;
		for(int k = 0; k <= 2 - n; k++)
		  for(int l = -k; l <= k; l++)
		    c_local[indice] += c_local_Ci[(n + k) * (n + k) + (n + k) + (m + l)] * Gnm[k * k + k + l];
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


  void L2P_low(Cell * Ci)
  {
    real_t r_local_Ci[9];
    complex_t c_local_Ci[9];

    for(int indice = 0; indice < 9; indice++)
      r_local_Ci[indice] = Ci->L[indice];

    real_2_complex(r_local_Ci, c_local_Ci, 2);

    for(Body * B = Ci->BODY; B != Ci->BODY + Ci->NBODY; B++)
      {
	complex_t Phi[4];

	for(int d = 0; d < 3; d++)
	  dX[d] = Ci->X[d] - B->X[d];

	complex_t Gnm[9];
	make_Gnm(&dX[0], &Gnm[0], 2);

	for(int i = 0; i < 9; i++)
	  Gnm[i] = std::conj(Gnm[i]);

	for(int n = 0; n <= 1; n++)
	  {
	    for(int m = 0; m <= n; m++)
	      {
		int indice = n * n + n + m;
		for(int k = 0; k <= 2 - n; k++)
		  for(int l = -k; l <= k; l++)
		    Phi[indice] += c_local_Ci[(n + k) * (n + k) + (n + k) + (m + l)] * Gnm[k * k + k + l];
	      }
	  }

	B->acc_old += sqrt(std::real(Phi[3]) * std::real(Phi[3]) +
			   std::imag(Phi[3]) * std::imag(Phi[3]) + std::real(Phi[2]) * std::real(Phi[2]));
      }
  }



}
#endif
