#ifndef kernel_best_lean_h
#define kernel_best_lean_h

#include "exafmm.h"
#include <iostream>
#include <math.h>		/*isnan, sqrt, atan2 */
#include <stdio.h>		/*printf */

#include "vectorclass/vectorclass.h"
#define NSIMD 8
#include "kernel_rotate.h"

#define M_SQRT2 1.41421356237309504880168872420969808
#define TINY 1.52587890625e-05

double dt_param = 0.01;

namespace exafmm
{

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

  void real_2_complex(real_t * real_arr, complex_t * complex_arr)
  {
    for(int n = 0; n <= P; n++)
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

  void complex_2_real(complex_t * complex_arr, real_t * real_arr)
  {
    for(int n = 0; n <= P; n++)
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
	int nii = ceil((float) ni / (float) NSIMD) * NSIMD;
	double Xi[nii] __attribute__ ((aligned(32)));
	double Yi[nii] __attribute__ ((aligned(32)));
	double Zi[nii] __attribute__ ((aligned(32)));

	for(int k = 0; k < ni; k++)
	  {
	    Xi[k] = Bi[k].X[0];
	    Yi[k] = Bi[k].X[1];
	    Zi[k] = Bi[k].X[2];
	  }

	Vec8d xi, yi, zi, r, mj;
	Vec8d invR, factor1, fac1, dx, dy, dz, r2;
	Vec8d ax, ay, az, pot;

	for(int i = 0; i < nii; i = i + NSIMD)
	  {
	    xi.load(Xi + i);
	    yi.load(Yi + i);
	    zi.load(Zi + i);

	    ax = 0, ay = 0, az = 0, pot = 0;

	    for(int j = 0; j < nj; j++)
	      {
		mj = Bj[j].q;
		dx = xi - Bj[j].X[0];
		dy = yi - Bj[j].X[1];
		dz = zi - Bj[j].X[2];

		r2 = dx * dx + dy * dy + dz * dz;

		mj = select(r2 > 0, mj, 0);
		r2 = select(r2 > 0, r2, 1e38);

		r = sqrt(r2);
		invR = 1 / r;
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
		    Bi[i + k].p += (real_t) pot[k];
#pragma omp atomic
		    Bi[i + k].F[0] -= (real_t) ax[k];
#pragma omp atomic
		    Bi[i + k].F[1] -= (real_t) ay[k];
#pragma omp atomic
		    Bi[i + k].F[2] -= (real_t) az[k];
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
		    dX[d] = Bi[i].X[d] - Bj[j].X[d];
		    dV[d] = Bi[i].V[d] - Bj[j].V[d];
		  }

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

#if SIMD_P2P
	  int nii = ceil((float) ni / (float) NSIMD) * NSIMD;

	  double Mi[nii] __attribute__ ((aligned(32)));
	  double Xi[nii] __attribute__ ((aligned(32)));
	  double Yi[nii] __attribute__ ((aligned(32)));
	  double Zi[nii] __attribute__ ((aligned(32)));
	  double VXi[nii] __attribute__ ((aligned(32)));
	  double VYi[nii] __attribute__ ((aligned(32)));
	  double VZi[nii] __attribute__ ((aligned(32)));

	  for(int k = 0; k < ni; k++)
	    {
	      Mi[k] = Bi[k].q;
	      Xi[k] = Bi[k].X[0];
	      Yi[k] = Bi[k].X[1];
	      Zi[k] = Bi[k].X[2];
	      VXi[k] = Bi[k].V[0];
	      VYi[k] = Bi[k].V[1];
	      VZi[k] = Bi[k].V[2];
	    }

	  Vec8d mi, xi, yi, zi, vxi, vyi, vzi, r2, v2, vdotdr2, r, mj;
	  Vec8d invR, factor1, fac1, dx, dy, dz, dvx, dvy, dvz;
	  Vec8d ax, ay, az, pot, timestep, tau, dtau;

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
		  dx = xi - Bj[j].X[0];
		  dy = yi - Bj[j].X[1];
		  dz = zi - Bj[j].X[2];
		  dvx = vxi - Bj[j].V[0];
		  dvy = vyi - Bj[j].V[1];
		  dvz = vzi - Bj[j].V[2];

		  r2 = dx * dx + dy * dy + dz * dz;
		  v2 = dvx * dvx + dvy * dvy + dvz * dvz + 1e-20;

		  mj = select(r2 > 0, mj, 0);
		  mi += mj;
		  r2 = select(r2 > 0, r2, 1e38);

		  vdotdr2 = (dx * dvx + dy * dvy + dz * dvz) / r2;

		  r = sqrt(r2);
		  invR = 1 / r;

		  tau = dt_param / M_SQRT2 * r * sqrt(r / mi);
		  dtau = 1.5 * tau * vdotdr2;
		  dtau = select(dtau < 1, dtau, 1);
		  tau /= (1 - dtau * 0.5);
		  timestep = min(tau, timestep);

		  tau = dt_param * r / sqrt(v2);
		  dtau = tau * vdotdr2 * (1 + mi / (v2 * r));
		  dtau = select(dtau < 1, dtau, 1);
		  tau /= (1 - dtau * 0.5);
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
		      Bi[i + k].p += (real_t) pot[k];
#pragma omp atomic
		      Bi[i + k].F[0] -= (real_t) ax[k];
#pragma omp atomic
		      Bi[i + k].F[1] -= (real_t) ay[k];
#pragma omp atomic
		      Bi[i + k].F[2] -= (real_t) az[k];

		      if(Bi[i + k].timestep > timestep[k])
			Bi[i + k].timestep = timestep[k];
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
		      dX[d] = Bi[i].X[d] - Bj[j].X[d];
		      dV[d] = Bi[i].V[d] - Bj[j].V[d];
		    }

		  real_t R2 = norm(dX);

		  if(R2 > 0)
		    {
		      real_t R = sqrt(R2);
		      real_t vdotdr2;
		      real_t v2 = norm(dV);

		      vdotdr2 = (dX[0] * dV[0] + dX[1] * dV[1] + dX[2] * dV[2]) / R2;

		      real_t tau = dt_param / M_SQRT2 * sqrt(R * R2 / (Bi[i].q + Bj[j].q));
		      real_t dtau = 1.5 * tau * vdotdr2;
		      if(dtau > 1.)
			dtau = 1.;
		      tau /= (1 - dtau / 2);
		      if(tau < timestep)
			timestep = tau;

		      if(v2 > 0)
			{
			  tau = dt_param * R / sqrt(v2);
			  dtau = tau * vdotdr2 * (1 + (Bi[i].q + Bj[j].q) / (v2 * R));
			  if(dtau > 1.)
			    dtau = 1.;
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
		  Bi[i].F[0] -= ax;
#pragma omp atomic
		  Bi[i].F[1] -= ay;
#pragma omp atomic
		  Bi[i].F[2] -= az;

		  if(Bi[i].timestep > timestep)
		    Bi[i].timestep = timestep;
		}

	    }
#endif
	}
    }

    void P2M(Cell * C)
    {

      real_t min_acc = 1e30;

      if(1)
	{
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
	}

      real_t max_r2 = 1e-20;
      complex_t c_multipole[NTERM];

      for(Body * B = C->BODY; B != C->BODY + C->NBODY; B++)
	{
	  for(int d = 0; d < 3; d++)
	    {
	      dX[d] = B->X[d] - C->X[d];
	    }

	  max_r2 = fmax(max_r2, norm(dX));

	  complex_t Gnm[NTERM];

	  make_Gnm(&dX[0], &Gnm[0], P);

	  for(int indice = 0; indice < NTERM; indice++)
	    c_multipole[indice] += (B->q * B->issource) * Gnm[indice];
	}

      real_t r_multipole[NTERM];
      complex_2_real(&c_multipole[0], &r_multipole[0]);

      for(int indice = 0; indice < NTERM; indice++)
	C->M[indice] += r_multipole[indice];

      C->min_acc = min_acc;
      C->R = sqrt(max_r2);

    }

    void M2M(Cell * Ci)
    {

      real_t min_acc = 1e30;

      if(1)
	{
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

	      Ci->min_acc = min_acc;
	      Ci->R = max_r;
	    }
	}

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

	  real_2_complex(r_multipole_Cj, c_multipole_Cj);

	  complex_t Gnm[NTERM];
	  make_Gnm(&dX[0], &Gnm[0], P);

	  for(int n = 0; n <= P; n++)
	    {
	      for(int m = -n; m <= n; m++)
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
      complex_2_real(c_multipole, r_multipole);

      for(int indice = 0; indice < NTERM; indice++)
	Ci->M[indice] += r_multipole[indice];

    }

    void L2L(Cell * Ci)
    {
      real_t r_local_Ci[NTERM];
      complex_t c_local_Ci[NTERM];

      for(int indice = 0; indice < NTERM; indice++)
	r_local_Ci[indice] = Ci->L[indice];

      real_2_complex(r_local_Ci, c_local_Ci);

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
	      for(int m = -n; m <= n; m++)
		{
		  int indice = n * n + n + m;
		  for(int k = 0; k <= P - n; k++)
		    for(int l = -k; l <= k; l++)
		      c_local[indice] +=
			c_local_Ci[(n + k) * (n + k) + (n + k) + (m + l)] * Gnm[k * k + k + l];
		}
	    }

	  real_t r_local[NTERM];

	  complex_2_real(c_local, r_local);

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

      real_2_complex(r_local_Ci, c_local_Ci);

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
	      for(int m = -n; m <= n; m++)
		{
		  int indice = n * n + n + m;
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
	  for(int n = 0; n <= P; n++)
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
	  for(int n = 0; n <= P; n++)
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
	  for(int n = 0; n <= P; n++)
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
	  for(int n = 0; n <= P; n++)
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
	  // Store the results in real - values form              
	  for(int indice = 0; indice < NTERM; indice++)
	    Ci->L[indice] += r_local[indice];
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

	  for(int n = 0; n <= P; n++)
	    {
	      int index_start = n * n + n;
	      for(int m = 1; m <= n; m++)
		r_multipole_Cj[index_start - m] *= -1;
	    }

	  real_2_complex(r_multipole_Cj, c_Cj);

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

	  complex_2_real(c_local, r_local);

	  for(int indice = 0; indice < NTERM; indice++)
	    Ci->L[indice] += r_local[indice];
	}
    }
  }
#endif
