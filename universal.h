/* -*-C-*-

   Copyright (C) 2015 Massachusetts Institute of Technology

   universal.c is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or (at
   your option) any later version.

   universal.c is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with universal.c; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
   USA.

*/

/* written by Jack Wisdom in June and July, 2015
   David M. Hernandez helped test and improve it.
   Further edits by David M. Hernandez, summer 2020.
   If you use universal.c please cite Wisdom and Hernandez (2015).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _STATE_
#define _STATE_
typedef struct {
  double x, y, z, xd, yd, zd, xdd, ydd, zdd;
} Binary;
#endif

/* M_PI does not exist on all systems */
#define UNIVERSAL_PI           3.14159265358979323846

#define SUCCESS 1
#define FAILURE 0
#define MAXITER 1000
#define TOL     1e-15

static inline void add_cs2(double* p, double* csp, double inp){
  double y = inp - *csp;

  double t = *p + y;
  *csp = (t - *p) - y;
  *p = t;
}


int solve_universal_newton(double kc, double r0, double beta, double b, double eta, double zeta, double h,
			   double *X, double *S2, double *C2)
{
  double x,xp;
  double g1, g2, g3, arg;
  double s2, c2, xnew;
  int count;
  double cc;

  xnew = *X;
  x = 0.;

  count = 0;
  do {
    xp = x;
    x = xnew;
    arg = b*x/2.0;
    s2 = sin(arg);
    c2 = cos(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/beta;
    g3 = (x - g1)/beta;
    cc = eta*g1 + zeta*g2;
    xnew = (h + (x*cc - (eta*g2 + zeta*g3)))/(r0 + cc);
    if(count++ > MAXITER) {
      return(FAILURE);
    }
  }while((xnew != x) && (xnew != xp));
  //while(fabs((x-xnew)/xnew) > TOL);
  //while((xnew != x) && (xnew != xp));
  //GLOBALV = count;
  //printf("Kepler %0.16le %0.16le %0.16le %d\n",xnew,x,xp,count);
  //printf("count=%d\n",count);
 
  /* compute the output */
  x = xnew;
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);
  
  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

int solve_universal_laguerre(double kc, double r0, double beta, double b, double eta, double zeta, double h, 
			     double *X, double *S2, double *C2)
{
  double x,xp=0.;
  double g1, g2, g3;
  double arg, s2, c2, xnew, dx;
  double f, fp, fpp, g0;
  int count;
  double cc;

  xnew = *X;
  x = 0.;

  count = 0;

  do {
    double c5=5.0, c16 = 16.0, c20 = 20.0;

    x = xnew;
    arg = b*x/2.0;
    s2 = sin(arg);
    c2 = cos(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/beta;
    g3 = (x - g1)/beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    fp = r0 + eta*g1 + zeta*g2;
    g0 = 1.0 - beta*g2;
    fpp = eta*g0 + zeta*g1;
    dx = -c5*f/(fp + sqrt(fabs(c16*fp*fp - c20*f*fpp)));
    xnew = x + dx;
    if(count++ > MAXITER) {
      printf("failed Laguerre\n");
      //exit(0);
      return(FAILURE);
    }
  } while(fabs(dx) > TOL*fabs(xnew));

  /* refine with a last newton */
  /*
    count = 0;
    printf("entering %0.15le %0.15le %0.15le\n",xnew,x,xp);
    printf("entering Laguerre Newton\n");
    do {
  */
  xp = x;
  x = xnew;
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);
  g1 = 2.0*s2*c2/b;
  g2 = 2.0*s2*s2/beta;
  g3 = (x - g1)/beta;
  cc = eta*g1 + zeta*g2;
  xnew = (h + (x*cc - (eta*g2 + zeta*g3)))/(r0 + cc);
    
  //if(count++ > MAXITER) {
  //printf("check second Newt %0.15le %0.15le %0.15le\n",xnew,x,xp);
  //exit(0);
  /*  printf("failed second time\n");
      exit(0);
      return(FAILURE);
      }
      } while((xnew != x) && (xnew != xp));
  */
  
  x = xnew;
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);

  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

int solve_universal_bisection(double kc, double r0, double beta, double b, double eta, double zeta, double h, 
			      double *X, double *S2, double *C2)
{
  double x=0.,xp=0.;
  double g1, g2, g3;
  double arg, s2, c2;
  double f;
  double X_min, X_max, X_per_period, invperiod;
  int count;
  double xnew;
  double err;

  xnew = *X;

  err = TOL*fabs(xnew);

  /* bisection limits due to Rein et al. */
  invperiod = b*beta/(2.*UNIVERSAL_PI*kc);
  X_per_period = 2.*UNIVERSAL_PI/b;
  X_min = X_per_period * floor(h*invperiod);
  X_max = X_min + X_per_period;
  xnew = (X_max + X_min)/2.;

  count = 0;
  do {
    xp = x;
    x = xnew;
    arg = b*x/2.0;
    s2 = sin(arg);
    c2 = cos(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/beta;
    g3 = (x - g1)/beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    if (f>=0.){
      X_max = x;
    }else{
      X_min = x;
    }
    xnew = (X_max + X_min)/2.;
    if(count++ > MAXITER) return(FAILURE);

  } while((xnew != x) && (xnew != xp));

  x = xnew;
  
  arg = b*x/2.0;
  s2 = sin(arg);
  c2 = cos(arg);

  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

double sign(double x)
{
  if(x > 0.0) {
    return(1.0);
  } else if(x < 0.0) {
    return(-1.0);
  } else {
    return(0.0);
  }
}

double cubic1(double a, double b, double c)
{
  double Q, R;
  double theta, A, B, x1;
  
  Q = (a*a - 3.0*b)/9.0;
  R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  if(R*R < Q*Q*Q) {
    theta = acos(R/sqrt(Q*Q*Q));

    printf("three cubic roots %.16le %.16le %.16le\n", 
	   -2.0*sqrt(Q)*cos(theta/3.0) - a/3.0,
	   -2.0*sqrt(Q)*cos((theta + 2.0*UNIVERSAL_PI)/3.0) - a/3.0,
	   -2.0*sqrt(Q)*cos((theta - 2.0*UNIVERSAL_PI)/3.0) - a/3.0);
    exit(-1);

    return(x1);

  } else {
    A = -sign(R)*pow(fabs(R) + sqrt(R*R - Q*Q*Q), 1./3.);
    if(A == 0.0) {
      B = 0.0;
    } else {
      B = Q/A;
    }
    x1 = (A + B) - a/3.0;

    return(x1);
  }
}


int solve_universal_parabolic(double kc, double r0, double beta, double b, double eta, double zeta, double h, 
			      double *X, double *S2, double *C2)
{
  double x;
  double s2, c2;

  b = 0.0;

  s2 = 0.0;
  c2 = 1.0;
  /*
    g1 = x;
    g2 = x*x/2.0;
    g3 = x*x*x/6.0;
    g = eta*g1 + zeta*g2;
  */
    
  /*   f = r0*x + eta*g2 + zeta*g3 - h; */
  /*   this is just a cubic equation */
  x = cubic1(3.0*eta/zeta, 6.0*r0/zeta, -6.0*h/zeta);

  s2 = 0.0;
  c2 = 1.0;
  /* g1 = 2.0*s2*c2/b; */
  /* g2 = 2.0*s2*s2/beta; */

  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

int solve_universal_hyperbolic_newton(double kc, double r0, double minus_beta, double b, double eta, double zeta, double h, 
				      double *X, double *S2, double *C2)
{
  double x;
  double g1, g2, g3;
  double arg, s2, c2, xnew;
  int count;
  double g;

  xnew = *X;

  count = 0;
  do {
    x = xnew;
    arg = b*x/2.0;
    if(fabs(arg)>200.0) return(FAILURE);
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    g = eta*g1 + zeta*g2;
    xnew = (x*g - eta*g2 - zeta*g3 + h)/(r0 + g);

    if(count++ > MAXITER) return(FAILURE);

  } while(fabs(x - xnew) > TOL*fabs(xnew));
  
  x = xnew;
  arg = b*x/2.0;
  s2 = sinh(arg);
  c2 = cosh(arg);

  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

int solve_universal_hyperbolic_laguerre(double kc, double r0, double minus_beta, double b, double eta, double zeta, double h, 
					double *X, double *S2, double *C2)
{
  double x;
  double g0, g1, g2, g3;
  double arg, s2, c2, xnew;
  int count;
  double f;

  xnew = *X;

  count = 0;
  do {
    double c5=5.0, c16 = 16.0, c20 = 20.0;
    double fp, fpp, dx;
    double den;

    x = xnew;
    arg = b*x/2.0;
    if(fabs(arg)>50.0) return(FAILURE);
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    fp = r0 + eta*g1 + zeta*g2;
    g0 = 1.0 + minus_beta*g2;
    fpp = eta*g0 + zeta*g1;
    den = (fp + sqrt(fabs(c16*fp*fp - c20*f*fpp)));
    if(den == 0.0) return(FAILURE);
    dx = -c5*f/den;
    xnew = x + dx;
    if(count++ > MAXITER) return(FAILURE);

  } while(fabs(x - xnew) > TOL*fabs(xnew));
  
  /* refine with a step of Newton */
  {
    double g;

    x = xnew;
    arg = b*x/2.0;
    if(fabs(arg)>200.0) return(FAILURE);
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    g = eta*g1 + zeta*g2;
    xnew = (x*g - eta*g2 - zeta*g3 + h)/(r0 + g);
  }

  x = xnew;
  arg = b*x/2.0;
  s2 = sinh(arg);
  c2 = cosh(arg);

  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

int solve_universal_hyperbolic_bisection(double kc, double r0, double minus_beta, double b, double eta, double zeta, double h, 
					 double *X, double *S2, double *C2)
{
  double x;
  double g1, g2, g3;
  double arg, s2, c2, xnew;
  int count;
  double X_min, X_max;
  double f;
  double err;

  xnew = *X;

  err = TOL*fabs(xnew);

  X_min = 0.5*xnew;
  X_max = 10.0*xnew;

  {
    double fmin, fmax;
    x = X_min;
    arg = b*x/2.0;
    if(fabs(arg)>200.0) return(FAILURE);
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    fmin = r0*x + eta*g2 + zeta*g3 - h;

    x = X_max;
    arg = b*x/2.0;
    if(fabs(arg)>200.0) {
      x = 200.0/(b/2.0);
      arg = 200.0;
    }
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    fmax = r0*x + eta*g2 + zeta*g3 - h;

    if(fmin*fmax > 0.0) {
      return(FAILURE);
    }
  }

  count = 0;
  do {
    x = xnew;
    arg = b*x/2.0;
    if(fabs(arg)>200.0) return(FAILURE);
    s2 = sinh(arg);
    c2 = cosh(arg);
    g1 = 2.0*s2*c2/b;
    g2 = 2.0*s2*s2/minus_beta;
    g3 = -(x - g1)/minus_beta;
    f = r0*x + eta*g2 + zeta*g3 - h;
    if (f>=0.){
      X_max = x;
    }else{
      X_min = x;
    }
    xnew = (X_max + X_min)/2.;
    if(count++ > MAXITER) return(FAILURE);

  } while(fabs(x - xnew) > err);
  
  x = xnew;
  arg = b*x/2.0;
  s2 = sinh(arg);
  c2 = cosh(arg);

  *X = x;
  *S2 = s2;
  *C2 = c2;

  return(SUCCESS);
}

double new_guess(double r0, double eta, double zeta, double dt)
{
  double s;
  double reta, disc;

  if(zeta != 0.0) {
    s = cubic1(3.0*eta/zeta, 6.0*r0/zeta, -6.0*dt/zeta);
  } else if(eta != 0.0) {
    reta = r0/eta;
    disc = reta*reta + 2.0*dt/eta;
    if(disc >= 0.0) {
      s = -reta + sqrt(disc);
    } else {
      s = dt/r0;
    }
  } else {
    s = dt/r0;
  }
  
  return(s);
}

int kepler_step_internal(double kc, double dt, double beta, double b,
			 Binary *s0, Binary *s,
			 double r0, double v2, double eta, double zeta)
{
  double r;
  double c2, s2;
  double G1, G2;
  double a, x;
  int flag;
  double c, ca, bsa;

  if(beta < 0.0) {
    double x0;

    x0 = new_guess(r0, eta, zeta, dt);

    x = x0;
    flag = solve_universal_hyperbolic_newton(kc, r0, -beta, b, eta, zeta, dt, &x, &s2, &c2);

    if(flag == FAILURE) {
      x = x0;
      //flag = solve_universal_hyperbolic_laguerre(kc, r0, -beta, b, eta, zeta, dt, &x, &s2, &c2);
      flag = solve_universal_hyperbolic_bisection(kc, r0, -beta, b, eta, zeta, dt, &x, &s2, &c2);
    }

    if(flag == FAILURE) {
      return(FAILURE);
    }

    a = kc/(-beta);
    G1 = 2.0*s2*c2/b;
    c = 2.0*s2*s2;
    G2 = c/(-beta);
    ca = c*a;
    r = r0 + eta*G1 + zeta*G2;
    /*
      double csr = 0;
      r = r0;
      add_cs2(&r, &csr, zeta*G2);
      double etap = eta*G1 - csr;
      r += etap;
    */
    bsa = (a/r)*(b/r0)*2.0*s2*c2;

  } else if(beta > 0.0) {
    double x0;
    double ff, fp;

    /* x0 = dt/r0; */
    double period = (2.*UNIVERSAL_PI*kc)/(b*beta);
    double ratio = dt / period;
    dt = (ratio - (int)(ratio)) * period;
      
    x0 = dt/r0;
    ff = zeta*x0*x0*x0 + 3.0*eta*x0*x0;
    fp = 3.0*zeta*x0*x0 + 6.0*eta*x0 + 6.0*r0;
    x0 -= ff/fp;

    x = x0;
    flag = solve_universal_newton(kc, r0, beta, b, eta, zeta, dt, &x, &s2, &c2);
    if(flag == FAILURE) {
      x = x0;
      //flag = solve_universal_laguerre(kc, r0, beta, b, eta, zeta, dt, &x, &s2, &c2);
      flag = solve_universal_bisection(kc, r0, beta, b, eta, zeta, dt, &x, &s2, &c2);
    }

    if(flag == FAILURE) {
      return(FAILURE);
    }

    a = kc/beta;
    G1 = 2.0*s2*c2/b;
    c = 2.0*s2*s2;
    G2 = c/beta;
    ca = c*a;
    r = r0 + eta*G1 + zeta*G2;
    /*
      double csr = 0;
      r = r0;
      add_cs2(&r, &csr, zeta*G2);
      double etap = eta*G1 - csr;
      r += etap;
    */
    bsa = (a/r)*(b/r0)*2.0*s2*c2;

  } else {
    
    x = dt/r0;

    flag = solve_universal_parabolic(kc, r0, beta, b, eta, zeta, dt, &x, &s2, &c2);
    if(flag == FAILURE) {
      exit(-1);
    }
    
    G1 = x;
    G2 = x*x/2.0;
    ca = kc*G2;
    r = r0 + eta*G1 + zeta*G2;
    /*
      double csr = 0;
      r = r0;
      add_cs2(&r, &csr, zeta*G2);
      double etap = eta*G1 - csr;
      r += etap;
    */
    bsa = kc*x/(r*r0);
  }
  
  double g, fdot, fhat, gdothat;

  /* f = 1.0 - (ca/r0); */
  fhat = -(ca/r0);
  g = eta*G2 + r0*G1;
  fdot = -bsa;
    
  /* gdot = 1.0 - (ca/r); */
  gdothat = -(ca/r);

  s->x = s0->x + (fhat*s0->x + g*s0->xd);
  s->y = s0->y + (fhat*s0->y + g*s0->yd);
  s->z = s0->z + (fhat*s0->z + g*s0->zd);
  s->xd = s0->xd + (fdot*s0->x + gdothat*s0->xd);
  s->yd = s0->yd + (fdot*s0->y + gdothat*s0->yd);
  s->zd = s0->zd + (fdot*s0->z + gdothat*s0->zd);

  /*
    s->x = (fhat*s0->x + g*s0->xd);
    s->y = (fhat*s0->y + g*s0->yd);
    s->z = (fhat*s0->z + g*s0->zd);
    s->xd = (fdot*s0->x + gdothat*s0->xd);
    s->yd = (fdot*s0->y + gdothat*s0->yd);
    s->zd = (fdot*s0->z + gdothat*s0->zd);*/

  return(SUCCESS);
    
}

void kepler_step_depth(double kc, double dt, double beta, double b,
                       Binary *s0, Binary *s, int depth,
		       double r0, double v2, double eta, double zeta)
{
  int flag;
  Binary ss;

  if(depth > 30) {
    printf("kepler depth exceeded\n");
    exit(-1);
  }

  flag = kepler_step_internal(kc, dt, beta, b, s0, s, r0, v2, eta, zeta);

  if(flag == FAILURE) {
    kepler_step_depth(kc, dt/4.0, beta, b, s0, &ss, depth+1, r0, v2, eta, zeta);

    r0 = sqrt(ss.x*ss.x + ss.y*ss.y + ss.z*ss.z);
    v2 = ss.xd*ss.xd + ss.yd*ss.yd + ss.zd*ss.zd;
    eta = ss.x*ss.xd + ss.y*ss.yd + ss.z*ss.zd;
    zeta = kc - beta*r0;

    kepler_step_depth(kc, dt/4.0, beta, b, &ss, s, depth+1, r0, v2, eta, zeta);

    r0 = sqrt(s->x*s->x + s->y*s->y + s->z*s->z);
    v2 = s->xd*s->xd + s->yd*s->yd + s->zd*s->zd;
    eta = s->x*s->xd + s->y*s->yd + s->z*s->zd;
    zeta = kc - beta*r0;

    kepler_step_depth(kc, dt/4.0, beta, b, s, &ss, depth+1, r0, v2, eta, zeta);

    r0 = sqrt(ss.x*ss.x + ss.y*ss.y + ss.z*ss.z);
    v2 = ss.xd*ss.xd + ss.yd*ss.yd + ss.zd*ss.zd;
    eta = ss.x*ss.xd + ss.y*ss.yd + ss.z*ss.zd;
    zeta = kc - beta*r0;

    kepler_step_depth(kc, dt/4.0, beta, b, &ss, s, depth+1, r0, v2, eta, zeta);
  }
}


void kepler_step(double kc, double dt, Binary *s0, Binary *s)
{
  double r0, v2, eta, beta, zeta, b;

  r0 = sqrt(s0->x*s0->x + s0->y*s0->y + s0->z*s0->z);
  v2 = s0->xd*s0->xd + s0->yd*s0->yd + s0->zd*s0->zd;
  eta = s0->x*s0->xd + s0->y*s0->yd + s0->z*s0->zd;
  beta = 2.0*kc/r0 - v2;
  zeta = kc - beta*r0;
  b = sqrt(fabs(beta));

  kepler_step_depth(kc, dt, beta, b, s0, s, 0, r0, v2, eta, zeta);
}




