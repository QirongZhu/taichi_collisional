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


typedef std::pair<unsigned int, unsigned int> int2;
typedef uint32_t HashType;

inline int CLZ(unsigned int z) { return __builtin_clz(z); }
inline int CLZ(unsigned long z) { return __builtin_clzl(z); }
inline int CLZ(unsigned long long z) { return __builtin_clzll(z); }
inline int sign(int x) { return (x > 0) - (x < 0); }

inline uint64_t expandBits(uint64_t v)
{
  v = (v * 0x000100000001u) & 0xFFFF00000000FFFFu;
  v = (v * 0x000000010001u) & 0x00FF0000FF0000FFu;
  v = (v * 0x000000000101u) & 0xF00F00F00F00F00Fu;
  v = (v * 0x000000000011u) & 0x30C30C30C30C30C3u;
  v = (v * 0x000000000005u) & 0x9249249249249249u;
  return v;
}

inline uint32_t expandBits(uint32_t v)
{
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

struct Node
{
  Node *left = nullptr;
  Node *right = nullptr;
  int idx;
  int index;
  int NBODY;
  int BODY;
  bool selected = false;
  bool isleaf;
  bool isLeaf() {return isleaf;}
};

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
    int delta(int x, int y);
    int2 determineRange(int idx);
    int findSplit(int2 range);
    void assignNode(int idx);
    void getBoundBox();
    void sortBodies();

    void setBodies(Bodies &bodies_);
    void buildRadixTree();

    void flagNode();
    void flagNode(Node *n, int index);

    void sumUpward(Node *Ci);
    void sumUpward();    
    void convertCells();

void upwardPass(Cell *Ci);
        //! Upward pass interface
        void upwardPass();

        void upwardPass_low(Cell *Ci);
        //! Upward pass interface
        void upwardPass_low();

        void M2L_rotate(Cell *Ci, Cell *Cj);

        // Set the tree to be either (0) KD tree (1) Recursive Coordinates Bisection tree
        void setType(int i);

        //! Recursive call to dual tree traversal for horizontal pass
        void horizontalPass(Cell *Ci, Cell *Cj);
        //! Horizontal pass interface
        void horizontalPass();

        void P2M(Cell *C);
        void M2M(Cell *Ci);
        void P2P(Cell *Ci, Cell *Cj);
        void L2L(Cell *Ci);
        void L2P(Cell *Ci);

        void P2M_low(Cell *C);
        void M2M_low(Cell *Ci);
        void M2L_low(Cell *Ci, Cell *Cj);
        void P2P_low(Cell *Ci, Cell *Cj);
        void L2L_low(Cell *Ci);

        void L2P_low(Cell *Ci);

        void printTree();


  private:
    Bodies bodies;

    std::vector<std::pair<HashType, unsigned int> > MortonIDs;
    //typedef std::vector<int2, tbb::detail::d1::cache_aligned_allocator<int2>> MortonIDs;
    typedef std::vector<Node, tbb::detail::d1::cache_aligned_allocator<Node>> Nodes;


    Nodes Leafs;
    Nodes Tree;

    Cells cells;

    Forces forces;

    std::pair<real_t, real_t> range_x;
    std::pair<real_t, real_t> range_y;
    std::pair<real_t, real_t> range_z;
    real_t R0;
    real_t X0[3];
    Node *root = nullptr;
    int start = 0;
  };

}
