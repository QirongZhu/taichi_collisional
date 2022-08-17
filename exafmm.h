#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace exafmm
{
  //! Basic type definitions
  typedef double real_t;	//!< Floating point type
  typedef std::complex < real_t > complex_t;	//!< Complex type

  //! Structure of bodies
  struct Body
  {
    real_t X[3];	     	//!< Position
    real_t q;			//!< Charge
    real_t F[3];		//!< Force
    real_t p;			//!< Potential
    real_t V[3];
    real_t acc_old;
    unsigned int Morton[3];     //!< Morton IDs for tree build
    bool operator<(const Body &rhs) const {
      unsigned int x_=(Morton[0]^rhs.Morton[0]);
      unsigned int y_=(Morton[1]^rhs.Morton[1]);
      unsigned int z_=(Morton[2]^rhs.Morton[2]);
      if((z_>x_||((z_^x_)<x_&&(z_^x_)<z_))&&(z_>y_||((z_^y_)<y_ &&(z_^y_)<z_)))
	return Morton[2]<rhs.Morton[2];
      if((y_>x_||((y_^x_)<x_ && (y_^x_)<y_)))
	return Morton[1]<rhs.Morton[1];
      return Morton[0]<rhs.Morton[0];
    }
    real_t timestep;
    real_t dtimestep;
    int index;
    bool issink;
    bool issource;
  };
  typedef std::vector < Body > Bodies;	//!< Vector of bodies

  //! Structure of cells
  struct Cell
  {
    real_t *M;
    real_t *L;
    real_t *Pn;
    int NCHILD;			//!< Number of child cells
    int NBODY;			//!< Number of descendant bodies
    Cell *CHILD;		//!< Pointer of first child cell
    Body *BODY;			//!< Pointer of first body
    real_t X[3];		//!< Cell center
    real_t R;			//!< Cell radius
    real_t min_acc;
#if EXAFMM_LAZY
    std::vector<Cell*>listM2L;	//!< M2L interaction list
    std::vector<Cell*>listP2P;	//!< P2P interaction list
#endif
    bool has_sink;
    bool has_source;
  };
  typedef std::vector < Cell > Cells;	//!< Vector of cells

  struct particle
  {
    real_t pos[3];
    real_t pos_e[3];
    real_t mass;
    real_t vel[3];
    real_t vel_e[3];
    real_t timestep;
    real_t acc[3];
    real_t pot;
    real_t postime;
    real_t acc_old;
    unsigned int id;
  };

  //! Global variables
  static char input_fname[200];
  unsigned int numBodies;
  int snapnum;

  real_t t_now;
  real_t force_accuracy = 1e-5;

  const real_t G = 1;
  //gravitaional constant with Msun, pc and km/s.

  struct sys
  {
    size_t n;
    struct particle *part;
    struct particle *last;
  };

  struct sys mainsys;

  double dt_param = 0.025;

  real_t PI2 = M_PI * M_PI;
  real_t PI4 = PI2 * PI2;
  real_t PI6 = PI4 * PI2;
  
  int P;			//!< Order of expansions
  int NTERM;			//!< Number of coefficients
  const int ncrit    = 192;		//!< Number of bodies per leaf cell
  const real_t theta = 0.65;		//!< Multipole acceptance criterion
  real_t dX[3], dV[3];		//!< Distance vector
#pragma omp threadprivate(dX, dV)	//!< Make global variables private
}
#endif
