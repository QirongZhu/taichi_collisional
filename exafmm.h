#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <omp.h>

#define ENDRUN(fmt, ...) {				\
    printf("ENDRUN at %s:%d ", __FILE__, __LINE__);	\
    printf(fmt, ## __VA_ARGS__);			\
    fflush(stdout);					\
    exit(-1);						\
  }

#include "./CGAL/Min_sphere_of_spheres_d.h"
class Ball {
private: // representation:
  double c[3]; // center in Eucliden coordinates
  double r;    // radius
public: // constructor:
  Ball() {}
  template<typename InputIterator>
  Ball(InputIterator from, double r) : r(r) {
    c[0] = *from;
    c[1] = *++from;
    c[2] = *++from;
  }
public: // accessors:
  double radius() const { return r; }
public: // iterator to iterate over the 3 coordinates:
  typedef const double *Coord_iterator;
  Coord_iterator begin_center() const { return c; }
};

struct Ball_traits {
  typedef Ball Sphere;
  static const int D=3;
  typedef double FT;
  typedef CGAL::Default_algorithm Algorithm;
  typedef CGAL::Tag_false Use_square_roots;
  typedef Sphere::Coord_iterator Cartesian_const_iterator;
  static Cartesian_const_iterator
  center_cartesian_begin(const Ball& b) {
    return b.begin_center();
  }
  static double radius(const Ball& b) {
    return b.radius();
  }
};
typedef CGAL::Min_sphere_of_spheres_d<Ball_traits> Minsphere;

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
    //real_t cell_mass;
#if EXAFMM_LAZY
    std::vector<Cell*>listM2L;	//!< M2L interaction list
    std::vector<Cell*>listP2P;	//!< P2P interaction list
#endif
    //    omp_lock_t *p2p_lock;
    //    omp_lock_t *m2l_lock;
    bool has_sink;
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
  real_t force_accuracy = 2.0e-06;

  real_t G = 1;			//0.004300710573170628; 
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
  int ncrit    = 250;		//!< Number of bodies per leaf cell
  real_t theta = 0.5;		//!< Multipole acceptance criterion
  real_t dX[3], dV[3];		//!< Distance vector
#pragma omp threadprivate(dX, dV)	//!< Make global variables private
}
#endif
