#include "common.h"

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
  bool isleaf;
  bool isLeaf() {return isleaf;}
};

namespace FMM
{
  class RadixTree
  {
  public:
    int delta(int x, int y);
    int2 determineRange(int idx);
    int findSplit(int2 range);
    void assignNode(int idx);
    void getBoundBox();
    void sortMontonID();

    void setBodies(Bodies &bodies_);
    void buildRadixTree();

    void traverse();

    void traverse(Node *n);

    void upwardPass(Node *Ci);
        //! Upward pass interface
    void upwardPass();

    void printTree();

  private:
    Bodies bodies;

    std::vector<std::pair<HashType, unsigned int> > sortedMortonCodes;
    typedef std::vector<Node> Nodes;

    //typedef std::vector<int2, tbb::detail::d1::cache_aligned_allocator<int2>> sortedMortonCodes;
    //typedef std::vector<Node, tbb::detail::d1::cache_aligned_allocator<Node>> Nodes;


    Nodes Leafs;
    Nodes Tree;

    std::pair<real_t, real_t> range_x;
    std::pair<real_t, real_t> range_y;
    std::pair<real_t, real_t> range_z;
    real_t R0;
    real_t X0[3];
    Node *root = nullptr;
    int start = 0;
  };

}