#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <complex>
#include <random>
#include <string>
#include <chrono>
#include <numeric>
#include <parallel/algorithm>
#include <memory>

#include <tbb/scalable_allocator.h>
#include <tbb/cache_aligned_allocator.h>

namespace FMM
{
    typedef double real_t;                  //!< Real type
    typedef std::complex<real_t> complex_t; //!< Complex type

    enum TreeType
    {
        kdtree,
        rcbtree,
        octree
    };

    struct Body
    {
        real_t X[3];
        real_t m;
        size_t index;

        int octant;

        Body(real_t x_ = 0, real_t y_ = 0, real_t z_ = 0, real_t m_ = 0, size_t id = 0)
        {
            X[0] = x_, X[1] = y_, X[2] = z_, m = m_, index = id;
        }
    };

    typedef std::vector<Body, tbb::detail::d1::cache_aligned_allocator<Body>> Bodies;
    // typedef std::vector<Body, std::allocator<Body>> Bodies;

    typedef Bodies::iterator Biter;

    struct Force
    {
        real_t F[3];
        real_t p;
        Force()
        {
            F[0] = 0, F[1] = 0, F[2] = 0, p = 0;
        }
    };
    typedef std::vector<Force, tbb::detail::d1::cache_aligned_allocator<Force>> Forces;
    //  typedef std::vector<Force, std::allocator<Force>> Forces;

    struct Cell
    {
        real_t X[3];
        real_t R;

        real_t *M;
        real_t *L;
        real_t *Pn;
        real_t min_acc;

        int CHILD;
        int NCHILD;
        size_t BODY;
        size_t NBODY;
        size_t index;

        bool isLeaf() { return CHILD == -1; }
    };

    typedef std::vector<Cell, tbb::detail::d1::cache_aligned_allocator<Cell>> Cells;
    // typedef std::vector<Cell, std::allocator<Cell>> Cells;

    const int P = EXPANSION;                             //!< Order of expansions
    const int NTERM = (EXPANSION + 1) * (EXPANSION + 1); //!< Number of coefficients
    const int ncrit = LEAFSIZE;                          //!< Number of bodies per leaf cell
    const real_t theta = 0.5;                            //!< Multipole acceptance criterion
    const int para_thres = 1500;

    static real_t dX[3];      //!< Distance vector
#pragma omp threadprivate(dX) //!< Make global variables private

    static real_t factorial_table[2 * EXPANSION + 8];

    static real_t factorial_coef[EXPANSION + 1][EXPANSION + 1];
    static real_t factorial_coef_inv[EXPANSION + 1][EXPANSION + 1];

    static real_t factorial_coef_oned[NTERM];
    static real_t factorial_coef_inv_oned[NTERM];

    static real_t combinator_coef[EXPANSION + 1][EXPANSION + 1];

}
