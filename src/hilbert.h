/* Hilbert encoding rotine taken from exafmm.

02/24/2023: Added a static array for faster bit operations; A new function to find ancestors

*/

#pragma once

#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <array>
#include <cmath>
#include <tuple>

#define EXAFMM_HILBERT 1 //! Set this to 0 for Morton
#define MAX_LEVEL 21 //! Max level of octree

namespace hilbert {

  static uint64_t levelOffsets[MAX_LEVEL + 1];

  //! Levelwise offset of Hilbert key
  inline uint64_t levelOffset(int level) {
    return (((uint64_t)1 << 3 * level) - 1) / 7;
  }

  void computeOffsets()
  {
    for (int i=0; i<=MAX_LEVEL; i++)
    {
        levelOffsets[i] = levelOffset(i);
    }
  }

  //! Get level from Hilbert key
  int getLevel(uint64_t i) {
    int level = -1;
    uint64_t offset = 0;
    while (i >= offset) {
      level++;
      offset += (uint64_t)1 << 3 * level;
    }
    return level;
  }

  //! Get parent's Hilbert key
  uint64_t getParent(uint64_t i) {
    int level = getLevel(i);
    return (i - levelOffset(level)) / 8 + levelOffset(level-1);
  }

  //! Get ancestor's key giving the level 
    uint64_t getAncestor(uint64_t i, int ans_level)
    {
        int level = MAX_LEVEL; 
        while (level > ans_level)
        {
            i = (i - levelOffsets[level]) / 8 + levelOffsets[level - 1];
            level--;
        }
        return i;
    }

  //! Get first child's Hilbert key
  uint64_t getChild(uint64_t i) {
    int level = getLevel(i);
    return (i - levelOffset(level)) * 8 + levelOffset(level+1);
  }

  //! Determine which octant the key belongs to
  int getOctant(uint64_t key, bool offset=true) {
    int level = getLevel(key);
    if (offset) key -= levelOffset(level);
    return key & 7;
  }

  //! Get Hilbert key from 3-D index
  uint64_t getKey(std::vector<int>& iX, int level, bool offset=true) {
#if EXAFMM_HILBERT
    int M = 1 << (level - 1);
    for (int Q=M; Q>1; Q>>=1) {
      int R = Q - 1;
      for (int d=0; d<3; d++) {
        if (iX[d] & Q) iX[0] ^= R;
        else {
          int t = (iX[0] ^ iX[d]) & R;
          iX[0] ^= t;
          iX[d] ^= t;
        }
      }
    }
    for (int d=1; d<3; d++) iX[d] ^= iX[d-1];
    int t = 0;
    for (int Q=M; Q>1; Q>>=1)
      if (iX[2] & Q) t ^= Q - 1;
    for (int d=0; d<3; d++) iX[d] ^= t;
#endif
    uint64_t i = 0;
    for (int l=0; l<level; l++) {
      i |= (iX[2] & (uint64_t)1 << l) << 2*l;
      i |= (iX[1] & (uint64_t)1 << l) << (2*l + 1);
      i |= (iX[0] & (uint64_t)1 << l) << (2*l + 2);
    }
    if (offset) i += levelOffset(level);
    return i;
  }

  //! Get 3-D index from Hilbert key
  std::tuple<int, int, int> get3DIndex(uint64_t i) {
    int level = getLevel(i);
    i -= levelOffset(level);
    int iX[3] = {0, 0, 0};
    for (int l=0; l<level; l++) {
      iX[2] |= (i & (uint64_t)1 << 3*l) >> 2*l;
      iX[1] |= (i & (uint64_t)1 << (3*l + 1)) >> (2*l + 1);
      iX[0] |= (i & (uint64_t)1 << (3*l + 2)) >> (2*l + 2);
    }
#if EXAFMM_HILBERT
    int N = 2 << (level - 1);
    int t = iX[2] >> 1;
    for (int d=2; d>0; d--) iX[d] ^= iX[d-1];
    iX[0] ^= t;
    for (int Q=2; Q!=N; Q<<=1) {
      int R = Q - 1;
      for (int d=2; d>=0; d--) {
        if (iX[d] & Q) iX[0] ^= R;
        else {
          t = (iX[0] ^ iX[d]) & R;
          iX[0] ^= t;
          iX[d] ^= t;
        }
      }
    }
#endif
    return std::make_tuple(iX[0], iX[1], iX[2]);
  }

  //! Get 3-D index from Hilbert key without level offset
  std::tuple<int, int, int> get3DIndex(uint64_t i, int level) {
    int iX[3] = {0, 0, 0};
    for (int l=0; l<level; l++) {
      iX[2] |= (i & (uint64_t)1 << 3*l) >> 2*l;
      iX[1] |= (i & (uint64_t)1 << (3*l + 1)) >> (2*l + 1);
      iX[0] |= (i & (uint64_t)1 << (3*l + 2)) >> (2*l + 2);
    }
#if EXAFMM_HILBERT
    int N = 2 << (level - 1);
    int t = iX[2] >> 1;
    for (int d=2; d>0; d--) iX[d] ^= iX[d-1];
    iX[0] ^= t;
    for (int Q=2; Q!=N; Q<<=1) {
      int R = Q - 1;
      for (int d=2; d>=0; d--) {
        if (iX[d] & Q) iX[0] ^= R;
        else {
          t = (iX[0] ^ iX[d]) & R;
          iX[0] ^= t;
          iX[d] ^= t;
        }
      }
    }
#endif
    return std::make_tuple(iX[0], iX[1], iX[2]);
  }

  //! Get 3-D index from coordinates
  std::array<int, 3> get3DIndex(std::vector<double>& X, int level) {

    double dx = 1.0 / (1LL << level);
    int iX[3];
    for (int d=0; d<3; d++) {
      iX[d] = static_cast<int> (std::floor((X[d] - 0.0 ) / dx));
    }
    return std::array<int, 3> {iX[0], iX[1], iX[2]};
  }

    //! Get 3-D index from coordinates
  std::array<int, 3> get3DIndex(std::array<double, 3>& X, int level) {

    double dx = 1.0 / (1LL << level);
    int iX[3];
    for (int d=0; d<3; d++) {
      iX[d] = static_cast<int> (std::floor((X[d] - 0.0 ) / dx));
    }
    return std::array<int, 3> {iX[0], iX[1], iX[2]};
  }

  //! Get coordinates from 3-D index
  std::array<double, 3> getCoordinates(int iX, int iY, int iZ, int level) {

    double dx = 1.0 / (1LL << level);

    double X[3];

    X[0] = (iX + 0.5) * dx;
    X[1] = (iX + 0.5) * dx;
    X[2] = (iX + 0.5) * dx;

    return std::array<double, 3> {X[0], X[1], X[2]};
  }
}
