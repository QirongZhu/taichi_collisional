#include "common.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

  typedef CGAL::Simple_cartesian<double> Kp;
  typedef CGAL::Min_sphere_of_points_d_traits_3<Kp, double> Traitsp;
  typedef CGAL::Min_sphere_of_spheres_d<Traitsp> Min_spherep;
  typedef Kp::Point_3 Pointp;

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

  typedef double FT;
  typedef CGAL::Cartesian<FT> K;
  typedef CGAL::Min_sphere_of_spheres_d_traits_3<K, FT> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
  typedef K::Point_3 Point;
  typedef Traits::Sphere Sphere;

namespace FMM
{

    real_t norm(real_t *X);
    int index(int n, int m);

    void real_2_complex(real_t *real_arr, complex_t *complex_arr, int order);
    void complex_2_real(complex_t *complex_arr, real_t *real_arr, int order);
    void make_Tnm(real_t *dX, complex_t *Tnm, int Order);
    void make_Gnm(real_t *dX, complex_t *Gnm, int Order);
    void make_Gnm_real(real_t *dX, real_t *Gnm, int Order);


    class Tree
    {
    public:
        static void initKernel();
        int findSplitAxis(Biter begin, Biter end, double &split_pos);
        Biter splitBodiesAlong(Biter begin, Biter end, unsigned int i, double split_pos);
        void sortBodiesAlong(Biter begin, Biter end, unsigned int i);
        void buildCells(Biter begin, Biter end, size_t index, int level);
        void allocateMultipoles();
        void buildTree();
        //void setBodies(size_t num);
        void setBodies(Bodies& bodies_);
        void preorder(size_t index);
        void preOrder();
        void upwardPass(Cell *Ci);
        //! Upward pass interface
        void upwardPass();
        // Set the tree to be either (0) KD tree (1) Recursive Coordinates Bisection tree
        void setType(int i);
        void P2M(Cell *C);
        void M2M(Cell *Ci);
        void P2P(Cell *Ci, Cell *Cj);
        void M2L_rotate(Cell *Ci, Cell *Cj);
        void L2L(Cell *Ci);
        void L2P(Cell *Ci);
        void printTree();

    private:
        Bodies bodies;
        Cells cells;

        std::vector<std::vector<real_t>> multipoles;
        std::vector<std::vector<real_t>> locals;
        std::vector<std::vector<real_t>> pns;

        TreeType treetype = kdtree;
        int dim = 3;
        int count = 0;
    };
}
